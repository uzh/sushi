SushiFabric::Application.configure do
  # Settings specified here will take precedence over those in config/application.rb

  # Code is not reloaded between requests
  config.cache_classes = true

  # Full error reports are disabled and caching is turned on
  config.consider_all_requests_local       = false
  config.action_controller.perform_caching = true

  # Disable Rails's static asset server (Apache or nginx will already do this)
  #config.serve_static_assets = false
  #config.serve_static_assets = true
  config.serve_static_files = true

  # Compress JavaScripts and CSS
  config.assets.compress = true

  # Don't fallback to assets pipeline if a precompiled asset is missed
  config.assets.compile = true

  # Generate digests for assets URLs
  config.assets.digest = true

  # Defaults to nil and saved in location specified by config.assets.prefix
  # config.assets.manifest = YOUR_PATH

  # Specifies the header that your server uses for sending files
  # config.action_dispatch.x_sendfile_header = "X-Sendfile" # for apache
  # config.action_dispatch.x_sendfile_header = 'X-Accel-Redirect' # for nginx

  # Force all access to the app over SSL, use Strict-Transport-Security, and use secure cookies.
  # config.force_ssl = true

  # See everything in the log (default is :info)
  # config.log_level = :debug

  # Prepend all log lines with the following tags
  # config.log_tags = [ :subdomain, :uuid ]

  # Use a different logger for distributed setups
  # config.logger = ActiveSupport::TaggedLogging.new(SyslogLogger.new)

  # Use a different cache store in production
  # config.cache_store = :mem_cache_store

  # Enable serving of images, stylesheets, and JavaScripts from an asset server
  # config.action_controller.asset_host = "http://assets.example.com"

  # Precompile additional assets (application.js, application.css, and all non-JS/CSS are already added)
  # config.assets.precompile += %w( search.js )

  # Disable delivery errors, bad email addresses will be ignored
  # config.action_mailer.raise_delivery_errors = false

  # Enable threaded mode
  # config.threadsafe!

  # Enable locale fallbacks for I18n (makes lookups for any locale fall back to
  # the I18n.default_locale when a translation can not be found)
  config.i18n.fallbacks = true

  # Send deprecation notices to registered listeners
  config.active_support.deprecation = :notify

  # Log the query plan for queries taking more than this (works
  # with SQLite, MySQL, and PostgreSQL)
  # config.active_record.auto_explain_threshold_in_seconds = 0.5

  config.logger = Logger.new("log/production.log", 5, 10 * 1024 * 1024)
  config.logger.level = Logger::ERROR
  config.log_level = :info
  config.eager_load = true

  def config.fgcz?
    @fgcz ||= (`hostname`.chomp =~ /fgcz/)
  end

  # fgcz
  if config.fgcz?
    #config.workflow_manager = "druby://fgcz-s-032:40001" # development
    #config.workflow_manager = "druby://fgcz-s-032:50001" # production
    #config.workflow_manager = "druby://fgcz-s-032:70001" # demo
    config.workflow_manager = "druby://fgcz-h-031:40001" # debian10 production
    #config.workflow_manager = "druby://fgcz-h-030:40001" # debian10 production
    config.scratch_dir = "/scratch"
    #config.gstore_dir = File.join(Dir.pwd, 'public/gstore/projects')
    config.gstore_dir = "/srv/gstore/projects" # production
    #config.gstore_dir = "/srv/GT/analysis/course_sushi/public/gstore/projects" # demo
    config.sushi_app_dir = Dir.pwd
    config.module_source = "/usr/local/ngseq/etc/lmod_profile"
    config.course_mode = false
    #ENV['PATH'] = "/usr/local/ngseq/packages/Dev/Python/3.6.8/bin/:/usr/local/ngseq/opt/GxTx_Scripts_in_Python3/g-bin/:" + ENV['PATH']
  end

end
