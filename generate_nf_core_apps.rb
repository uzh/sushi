require_relative 'master/lib/nf_core_app_factory'

puts "Generating nf-core App files..."
NfCoreAppFactory.generate_app_files
puts "Done."
