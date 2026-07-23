require 'csv'
require 'digest'

# Core logic for the machine-callable registration API
# (app/controllers/api/v1/datasets_controller.rb).
#
# Receives the dataset manifest as CONTENT (never a server-side caller path),
# parses/validates it in memory, and performs a SUSHI-local, single-transaction
# insert of the data_set + samples. Idempotency is keyed on a composed
# fingerprint stored in data_sets.registration_key.
#
# Authn/authz (token validity + project scope) is enforced by the controller
# before these methods run; this service assumes the caller is already
# authorized for the given project_number.
class DatasetRegistrationService
  # --- fingerprint -------------------------------------------------------

  # Deterministic canonical form of the manifest text: normalize line endings,
  # right-strip each line, drop trailing blank lines. Row ORDER is preserved
  # (sample order is semantically significant).
  def self.canonical(tsv)
    lines = tsv.to_s.gsub(/\r\n?/, "\n").split("\n", -1).map(&:rstrip)
    lines.pop while lines.any? && lines.last.empty?
    lines.join("\n")
  end

  # Composed, caller-context idempotency key (API INV-4 / R-INV-2).
  # The same content under a different project/name/parent is a DIFFERENT
  # registration.
  def self.registration_key(tsv:, project_number:, name:, parent_id:)
    payload = [canonical(tsv), project_number.to_s, name.to_s, parent_id.to_s].join("\x00")
    Digest::SHA256.hexdigest(payload)
  end

  # --- validation --------------------------------------------------------

  # Parse + validate the manifest and pre-flight parent/idempotency without
  # mutating anything. Returns { ok:, checks:[], errors:[] }.
  def self.validate(dataset_tsv:, project_number:, name:, parent_id: nil, order_id: nil)
    checks = []
    errors = []

    rows = parse_manifest(dataset_tsv)
    if rows.nil?
      errors << { check: "manifest", detail: "cannot parse dataset_tsv as TSV" }
      checks << { check: "manifest", status: "failed" }
      return { ok: false, checks: checks, errors: errors }
    end
    if rows.length < 2
      errors << { check: "manifest", detail: "need one header line and at least one sample line" }
      checks << { check: "manifest", status: "failed" }
      return { ok: false, checks: checks, errors: errors }
    end
    headers = rows.first
    unless headers.any? { |h| h.to_s.include?("[File]") }
      errors << { check: "manifest", detail: "no [File] tag column present" }
      checks << { check: "manifest", status: "failed" }
      return { ok: false, checks: checks, errors: errors }
    end
    checks << { check: "manifest", status: "ok" }

    # parent integrity (INV-10): must exist and share the project.
    if parent_id.nil? || parent_id.to_s.empty?
      checks << { check: "parent", status: "not_applicable" }
    else
      parent = DataSet.find_by(id: parent_id)
      if parent.nil?
        checks << { check: "parent", status: "failed", detail: "parent_id #{parent_id} not found" }
        errors << { check: "parent", detail: "parent_id #{parent_id} not found" }
      elsif parent.project&.number.to_i != project_number.to_i
        checks << { check: "parent", status: "failed", detail: "parent_id #{parent_id} belongs to a different project" }
        errors << { check: "parent", detail: "parent_id #{parent_id} belongs to a different project" }
      else
        checks << { check: "parent", status: "ok" }
      end
    end

    # idempotency lookup (no mutation).
    key = registration_key(tsv: dataset_tsv, project_number: project_number, name: name, parent_id: parent_id)
    existing = DataSet.find_by(registration_key: key)
    if existing
      checks << { check: "idempotency", status: "already_registered", data_set_id: existing.id }
    else
      checks << { check: "idempotency", status: "new", data_set_id: nil }
    end

    { ok: errors.empty?, checks: checks, errors: errors }
  end

  # --- registration ------------------------------------------------------

  # Register the dataset. Returns one of:
  #   { http: 200, body: { data_set_id:, idempotent_replay: } }
  #   { http: 422, body: { error:, checks:, errors: } }
  #   { http: 409, body: { error: } }
  # owner_login (optional): the acting user's LDAP login (set by the controller
  # for a user-principal token). When present, the new data_set is owned by that
  # User (created on first use), so the UI "Who" shows the person rather than the
  # nil-user fallback. Static tokens pass nil (no person) and keep the fallback.
  def self.register(dataset_tsv:, project_number:, name:, parent_id: nil, order_id: nil, owner_login: nil)
    v = validate(dataset_tsv: dataset_tsv, project_number: project_number, name: name,
                 parent_id: parent_id, order_id: order_id)
    return { http: 422, body: { error: "invalid", checks: v[:checks], errors: v[:errors] } } unless v[:ok]

    key = registration_key(tsv: dataset_tsv, project_number: project_number, name: name, parent_id: parent_id)

    existing = DataSet.find_by(registration_key: key)
    return { http: 200, body: { data_set_id: existing.id, idempotent_replay: true } } if existing

    rows    = parse_manifest(dataset_tsv)
    headers = rows.shift

    begin
      data_set_id = nil
      ActiveRecord::Base.transaction do
        project = Project.find_by_number(project_number) || Project.create!(number: project_number.to_i)

        data_set = DataSet.new(name: name, project: project, registration_key: key)
        data_set.parent_id = parent_id if parent_id.present?
        data_set.user = User.find_or_create_by(login: owner_login) if owner_login.present?

        rows.each do |row|
          sample_hash = {}
          headers.each_with_index { |header, i| sample_hash[header] = row[i] }
          data_set.samples.build(key_value: sample_hash.to_s)
        end

        data_set.md5 = data_set.md5hexdigest
        data_set.save!
        data_set_id = data_set.id
      end
      { http: 200, body: { data_set_id: data_set_id, idempotent_replay: false } }
    rescue ActiveRecord::RecordNotUnique
      # Concurrent duplicate won the unique (scope,key) race → treat as replay.
      existing = DataSet.find_by(registration_key: key)
      if existing
        { http: 200, body: { data_set_id: existing.id, idempotent_replay: true } }
      else
        { http: 409, body: { error: "conflict" } }
      end
    end
  end

  # --- set-once bfabric id ----------------------------------------------

  # Returns { http:, body: }. Set-once: same value is idempotent (200),
  # a different value is a conflict (409).
  def self.set_bfabric_id(data_set, bfabric_id)
    bfabric_id = bfabric_id.to_i
    if data_set.bfabric_id.present?
      if data_set.bfabric_id.to_i == bfabric_id
        { http: 200, body: { ok: true } }
      else
        { http: 409, body: { error: "bfabric_id already set to #{data_set.bfabric_id}" } }
      end
    else
      data_set.update!(bfabric_id: bfabric_id)
      { http: 200, body: { ok: true } }
    end
  end

  # --- deregister (compensation) ----------------------------------------

  def self.deregister(data_set)
    data_set.samples.destroy_all
    data_set.destroy!
    { http: 200, body: { ok: true, state: "COMPENSATED" } }
  end

  # --- helpers -----------------------------------------------------------

  # Parse TSV content into rows; returns nil on parse failure.
  def self.parse_manifest(content)
    CSV.parse(content.to_s, col_sep: "\t")
  rescue CSV::MalformedCSVError, ArgumentError
    nil
  end
end
