# Be sure to restart your server when you modify this file.

# Your secret key for verifying the integrity of signed cookies.
# If you change this key, all old signed cookies will become invalid!
# Make sure the secret is at least 30 characters and all random,
# no regular words or you'll be exposed to dictionary attacks.
secret_token = ENV['SUSHI_SECRET_TOKEN']
if secret_token.nil? || secret_token.length < 30
  raise 'SUSHI_SECRET_TOKEN env var must be set to a random string of at least 30 characters. Generate with: ruby -rsecurerandom -e "puts SecureRandom.hex(64)"'
end
SushiFabric::Application.config.secret_token = secret_token
