namespace :api_token do
  # Interim admin flow to mint a registration-API token. The raw token is
  # printed ONCE; only its salted hash is stored. A full issue/rotate/revoke
  # UI is backlog.
  #
  # static (default): scope is a fixed list of project numbers.
  #   rake api_token:issue NAME=registrar-test SCOPE=2680,3186 TTL_DAYS=90
  # user: bound to an LDAP login; authorized live against that login's current
  # FGCZ project membership. login + a bounded TTL are mandatory; SCOPE ignored.
  #   rake api_token:issue PRINCIPAL=user NAME=hubert-reg LOGIN=hubert TTL_DAYS=90
  desc "Issue an API token (NAME=, PRINCIPAL=static|user, SCOPE=comma,projects (static), LOGIN= (user), TTL_DAYS=)"
  task issue: :environment do
    name      = ENV["NAME"] or abort("NAME is required")
    principal = (ENV["PRINCIPAL"] || "static").strip
    ttl       = ENV["TTL_DAYS"]

    begin
      if principal == "user"
        raw, record = ApiToken.issue(name: name, ttl_days: ttl,
                                     principal: "user", login: ENV["LOGIN"])
      else
        scope = ENV["SCOPE"].to_s.split(",").map(&:strip).reject(&:empty?).map(&:to_i)
        abort("SCOPE is required (comma-separated project numbers)") if scope.empty?
        raw, record = ApiToken.issue(name: name, scope: scope, ttl_days: ttl)
      end
    rescue ArgumentError => e
      abort("cannot issue token: #{e.message}")
    end

    detail = record.user? ? "login=#{record.login}" : "scope=#{record.scope.inspect}"
    puts "Issued API token id=#{record.id} name=#{record.name} principal=#{record.principal} " \
         "#{detail} expires_at=#{record.expires_at || 'never'}"
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
