class SubmitJob < ApplicationJob
  queue_as :default

  # Safely evaluate parameter values - only allow numeric/boolean/array literals
  def safe_eval(value)
    return value if value.nil?
    str = value.to_s.strip
    
    # Boolean
    return true if str == 'true'
    return false if str == 'false'
    
    # Integer
    return str.to_i if str =~ /\A-?\d+\z/
    
    # Float
    return str.to_f if str =~ /\A-?\d+\.\d+\z/
    
    # Array literal (simple cases only)
    if str =~ /\A\[.*\]\z/
      begin
        return JSON.parse(str)
      rescue
        return value
      end
    end
    
    # Default: return as string
    value
  end

  def perform(params)
    class_name = params[:class_name]
    # For dynamically defined nf-core apps, ensure they are registered in this process
    unless Object.const_defined?(class_name)
      begin
        require 'nf_core_app_factory'
        NfCoreAppFactory.register_dynamic_apps
      rescue LoadError
        # nf_core_app_factory not available, continue with normal require
      end
    end
    # Still try require if class is not defined (for static .rb files)
    require class_name unless Object.const_defined?(class_name)
    sushi_app = eval(class_name).new
    sushi_app.logger = logger
    sushi_app.sushi_server = eval(SushiFabric::Application.config.sushi_server_class).new
    sushi_app.user = params[:user]
    sushi_app.next_dataset_name = params[:next_dataset_name]
    sushi_app.next_dataset_comment = params[:next_dataset_comment]
    sushi_app.params['sushi_app'] = class_name
    params[:parameters].each do |key, value|
      sushi_app.params[key] = if sushi_app.params.data_type(key) == String
                                       value
                                     else
                                       safe_eval(value)
                                     end
    end
    sushi_app.project = params[:project]
    sushi_app.dataset_sushi_id = params[:data_set_id]
    sushi_app.current_user = params[:current_user]
    sushi_app.off_bfabric_registration = params[:off_bfabric_registration]
    
    # Save project defaults if requested
    if params[:save_as_default]
      temp_file = sushi_app.save_project_defaults
      if temp_file && File.exist?(temp_file)
        temp_dir = File.dirname(temp_file)
        begin
          # Copy the temporary file to the project folder using copy_commands
          project_dir = File.join(SushiFabric::GSTORE_DIR, params[:project])
          
          # Use copy_commands to copy the file to gstore (force to overwrite existing file)
          copy_cmds = sushi_app.copy_commands(temp_file, project_dir, 'force')
          copy_cmds.each do |command|
            logger.info "Copying project defaults: #{command}"
            unless system command
              logger.error "Failed to copy project defaults file from #{temp_file} to #{project_dir}"
            end
          end
          
          # Clean up temporary directory and file
          FileUtils.rm_rf(temp_dir) if temp_dir && File.exist?(temp_dir)
        rescue => e
          logger.error "Error saving project defaults: #{e.message}"
          # Clean up temporary directory even if there's an error
          FileUtils.rm_rf(temp_dir) if temp_dir && File.exist?(temp_dir)
        end
      end
    end
    
    if params[:submit_type] == "Submit"
      sushi_app.run
    elsif params[:submit_type] == "MockRun"
      sushi_app.mock_run
    end
  end
end
