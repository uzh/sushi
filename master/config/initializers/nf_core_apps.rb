# Load nf-core apps dynamically at Rails boot
require_relative '../../lib/nf_core_app_factory'
require_relative '../../lib/nf_core_info_fetcher'

Rails.application.config.after_initialize do
  begin
    NfCoreAppFactory.register_dynamic_apps
    Rails.logger.info "Registered nf-core apps: #{NfCoreAppFactory.registered_class_names.join(', ')}"
  rescue => e
    Rails.logger.error "Failed to register nf-core apps: #{e.message}"
  end
end
