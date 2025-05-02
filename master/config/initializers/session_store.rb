# Be sure to restart your server when you modify this file.

SushiFabric::Application.config.session_store :cookie_store,
  key: '_sushi-fabric_session',
  domain: :all,
  tld_length: 3,
  secure: false,  # false for development
  same_site: :lax  # Security setting

# Use the database for sessions instead of the cookie-based default,
# which shouldn't be used to store highly confidential information
# (create the session table with "rails generate session_migration")
# SushiFabric::Application.config.session_store :active_record_store
