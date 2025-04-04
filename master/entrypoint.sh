#!/bin/bash
set -e

# Remove a potentially pre-existing server.pid for Rails.
rm -f /app/tmp/pids/server.pid

# Install dependencies if needed
bundle check || bundle install

# Initialize database if it doesn't exist
if [ "$RAILS_ENV" = "development" ] && [ ! -f /app/db/development.sqlite3 ]; then
  bundle exec rails db:create
  bundle exec rails db:migrate
  bundle exec rails db:seed
elif [ "$RAILS_ENV" = "production" ] && [ ! -f /app/db/production.sqlite3 ]; then
  bundle exec rails db:create
  bundle exec rails db:migrate
  bundle exec rails db:seed
fi

# Make sure the gstore projects directory exists and is accessible
mkdir -p /srv/gstore/projects
#chmod 755 /srv/gstore/projects

# Then exec the container's main process (what's set as CMD in the Dockerfile).
exec "$@" 
