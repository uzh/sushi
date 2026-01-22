# Load nf-core apps dynamically at Rails boot
puts "=== nf_core_apps.rb initializer loading ==="
$stdout.flush

require_relative '../../lib/nf_core_info_fetcher'
require_relative '../../lib/nf_core_app_factory'

puts "=== nf_core_apps.rb: Required files loaded ==="
$stdout.flush

Rails.application.config.after_initialize do
  puts "=== nf_core_apps.rb: after_initialize block executing ==="
  $stdout.flush
  begin
    puts "NfCoreApps: Starting to register nf-core apps..."
    $stdout.flush
    NfCoreAppFactory.register_dynamic_apps
    puts "NfCoreApps: Successfully registered #{NfCoreAppFactory.registered_class_names.size} nf-core apps"
    $stdout.flush
  rescue => e
    puts "NfCoreApps: Failed to register nf-core apps: #{e.message}"
    puts e.backtrace.first(10).join("\n")
    $stdout.flush
  end
end
