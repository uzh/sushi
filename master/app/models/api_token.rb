require 'digest'
require 'securerandom'

# Per-caller bearer token for the machine-callable registration API
# (app/controllers/api/v1/datasets_controller.rb).
#
# Secrecy at rest: only a salted SHA-256 hash of the raw token is stored; the
# raw token is shown exactly once at issue time and never persisted. Each token
# carries an explicit project scope and an optional expiry, and is revocable.
class ApiToken < ActiveRecord::Base
  serialize :scope, Array

  # Salt the hash with the app's secret_key_base so a leaked api_tokens table
  # is not directly reversible without the server secret.
  def self.salt
    Rails.application.secret_key_base.to_s
  end

  def self.digest(raw)
    Digest::SHA256.hexdigest(salt + raw.to_s)
  end

  # Issue a new token. Returns [raw_token, record]. Persist nothing else; the
  # raw value cannot be recovered afterwards.
  def self.issue(name:, scope:, ttl_days: nil)
    raw = SecureRandom.urlsafe_base64(32)
    record = create!(
      name:       name,
      token_hash: digest(raw),
      scope:      Array(scope).map(&:to_i),
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

  def active?
    !revoked? && !expired?
  end

  def revoked?
    revoked_at.present?
  end

  def expired?
    expires_at.present? && expires_at <= Time.now
  end

  def in_scope?(project_number)
    Array(scope).map(&:to_i).include?(project_number.to_i)
  end
end
