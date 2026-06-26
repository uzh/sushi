namespace :api_token do
  # Interim admin flow to mint a registration-API token. The raw token is
  # printed ONCE; only its salted hash is stored. A full issue/rotate/revoke
  # UI is backlog.
  #
  #   rake api_token:issue NAME=registrar-test SCOPE=2680,3186 TTL_DAYS=90
  desc "Issue an API token (NAME=, SCOPE=comma,separated,project,numbers, TTL_DAYS=optional)"
  task issue: :environment do
    name  = ENV["NAME"] or abort("NAME is required")
    scope = ENV["SCOPE"].to_s.split(",").map(&:strip).reject(&:empty?).map(&:to_i)
    abort("SCOPE is required (comma-separated project numbers)") if scope.empty?
    ttl   = ENV["TTL_DAYS"]

    raw, record = ApiToken.issue(name: name, scope: scope, ttl_days: ttl)
    puts "Issued API token id=#{record.id} name=#{record.name} scope=#{scope.inspect} " \
         "expires_at=#{record.expires_at || 'never'}"
    puts "RAW TOKEN (shown once, store securely):"
    puts raw
  end

  desc "Revoke an API token by id (ID=)"
  task revoke: :environment do
    id = ENV["ID"] or abort("ID is required")
    token = ApiToken.find(id)
    token.update!(revoked_at: Time.now)
    puts "Revoked API token id=#{token.id} name=#{token.name}"
  end
end
