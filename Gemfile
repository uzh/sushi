source 'https://rubygems.org'

gem 'rails', '3.2.9'

# Bundle edge Rails instead:
# gem 'rails', :git => 'git://github.com/rails/rails.git'

gem 'sqlite3'


# Gems used only for assets and not required
# in production environments by default.
group :assets do
  gem 'sass-rails',   '~> 3.2.3'
  gem 'coffee-rails', '~> 3.2.1'

  # See https://github.com/sstephenson/execjs#readme for more supported runtimes
  # gem 'therubyracer', :platforms => :ruby

  gem 'uglifier', '>= 1.0.3'
end

gem 'jquery-rails'

# To use ActiveModel has_secure_password
# gem 'bcrypt-ruby', '~> 3.0.0'

# To use Jbuilder templates for JSON
# gem 'jbuilder'

# Use unicorn as the app server
# gem 'unicorn'

# Deploy with Capistrano
# gem 'capistrano'

# To use debugger
# gem 'debugger'
gem 'devise'
gem 'devise_ldap_authenticatable'
gem 'hpricot'
gem 'ruby_parser'

# Options for FGCZ gems:
# A) gem 'fgcz', :bzr => 'bzr+ssh://fgcz-s-034.uzh.ch/usr/local/ngseq/repo/Gems/fgcz'
#    but bzr source doesn't exist! (only git)
# B) gem 'fgcz', :path => '/srv/SushiFabric/Gems/fgcz'
# C) source 'https://gems-fgcz.uzh.ch:8888'
#    gem 'fgcz'
gem 'savon'
if `hostname`.chomp =~ /fgcz-s-034/
gem 'fgcz', :path => '/srv/SushiFabric/Gems/fgcz'
end

group :test do
  gem "rspec"
  gem "rspec-rails"
  gem 'simplecov'
  gem 'simplecov-rcov'
end

gem 'workflow_manager'
