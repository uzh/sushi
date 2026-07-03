require 'digest'
require 'securerandom'

# Per-caller bearer token for the machine-callable registration API
# (app/controllers/v1/datasets_controller.rb).
#
# Secrecy at rest: only a salted SHA-256 hash of the raw token is stored; the
# raw token is shown exactly once at issue time and never persisted.
#
# A token has exactly one principal (design v0.7):
#   - `static`: authorized against the stored project-number array (`scope`),
#     frozen at issue. Existing/legacy behavior; unchanged (P-INV-8).
#   - `user`: bound to a non-blank LDAP `login` and a mandatory bounded TTL;
#     authorized live, per request, against the login's current FGCZ project
#     membership (no cache, W=0). `scope` is unused.
class ApiToken < ActiveRecord::Base
  serialize :scope, Array

  # Mandatory upper bound on a `user` token's lifetime (design v0.7 P-INV-10).
  MAX_USER_TOKEN_TTL_DAYS = 90

  # Raised when the live project-membership resolver cannot answer (transport
  # failure). The controller maps this to 503, distinct from an authorization
  # denial (403). Fail-closed: the request never proceeds.
  class ResolverUnavailable < StandardError; end

  # Salt the hash with the app's secret_key_base so a leaked api_tokens table
  # is not directly reversible without the server secret.
  def self.salt
    Rails.application.secret_key_base.to_s
  end

  def self.digest(raw)
    Digest::SHA256.hexdigest(salt + raw.to_s)
  end

  # Issue a new token. Returns [raw_token, record]. The raw value cannot be
  # recovered afterwards.
  #
  # principal: "static" (default) requires a non-empty scope; "user" requires a
  # non-blank login and a TTL within MAX_USER_TOKEN_TTL_DAYS (scope is ignored).
  def self.issue(name:, scope: [], ttl_days: nil, principal: "static", login: nil)
    principal = principal.to_s
    case principal
    when "user"
      raise ArgumentError, "login is required for a user token" if login.to_s.strip.empty?
      ttl_raw = ttl_days.to_s.strip
      raise ArgumentError, "a user token requires TTL_DAYS (mandatory bounded TTL)" if ttl_raw.empty?
      # Strict integer: reject "90abc"/"90.9" rather than silently truncating.
      raise ArgumentError, "TTL_DAYS must be a positive integer" unless ttl_raw.match?(/\A\d+\z/)
      ttl_days = ttl_raw.to_i
      if ttl_days <= 0 || ttl_days > MAX_USER_TOKEN_TTL_DAYS
        raise ArgumentError, "TTL_DAYS must be between 1 and #{MAX_USER_TOKEN_TTL_DAYS} for a user token"
      end
      scope = []
    when "static"
      # A static token authorizes only its explicit scope; an empty scope would
      # authorize nothing, so require it (the `scope:` default exists for the
      # user branch, which ignores scope).
      if Array(scope).map { |x| x.to_s.strip }.reject(&:empty?).empty?
        raise ArgumentError, "scope is required for a static token"
      end
    else
      raise ArgumentError, "unknown principal #{principal.inspect} (expected static|user)"
    end

    raw = SecureRandom.urlsafe_base64(32)
    record = create!(
      name:       name,
      token_hash: digest(raw),
      scope:      Array(scope).map(&:to_i),
      principal:  principal,
      login:      (principal == "user" ? login.to_s.strip : nil),
      expires_at: ttl_days ? Time.now + ttl_days.to_i.days : nil
    )
    [raw, record]
  end

  # Look up an active token by its raw value. Fail-closed: returns nil for any
  # missing/unknown/expired/revoked token.
  def self.authenticate(raw)
    return nil if raw.to_s.empty?
    token = find_by(token_hash: digest(raw))
    return nil unless token && token.active?
    token
  end

  def user?
    principal.to_s == "user"
  end

  def static?
    !user?
  end

  def active?
    return false if revoked? || expired?
    # A user token's TTL is mandatory; a null expiry must never authenticate
    # (defense in depth — issuance already enforces this). Design v0.7 step 2.
    return false if user? && expires_at.nil?
    true
  end

  def revoked?
    revoked_at.present?
  end

  def expired?
    expires_at.present? && expires_at <= Time.now
  end

  # Static-principal membership test (unchanged). Not used for `user` tokens;
  # user membership is tested against the live-resolved set (see allowed_projects).
  def in_scope?(project_number)
    Array(scope).map(&:to_i).include?(project_number.to_i)
  end

  # The set of project numbers this token may currently act on.
  #   - static: the stored scope array.
  #   - user:   the login's current FGCZ project membership, resolved live
  #             (W=0). An inactive/unknown login yields the empty set (→ 403).
  #
  # Raises ResolverUnavailable when the resolver *call* fails (→ 503). NOTE:
  # FGCZ.get_user_projects2 swallows most LDAP failures into an empty array, so a
  # silent backend outage degrades to []→403 rather than 503; either way the
  # request is denied (fail-closed), only the HTTP granularity is lost. The
  # rescue is scoped to the resolver call alone so parsing bugs are NOT masked as
  # infrastructure errors.
  def allowed_projects
    return Array(scope).map(&:to_i) unless user?

    raw =
      begin
        FGCZ.get_user_projects2(login)
      rescue => e
        raise ResolverUnavailable, "#{e.class}: #{e.message}"
      end

    Array(raw).map { |p| p.to_s.sub(/\Ap/i, "").to_i }.select(&:positive?)
  end
end
