#!/usr/bin/env ruby
# encoding: utf-8
# Version = '20250620-102347'

require 'csv'
require 'fileutils'
require 'yaml'
require 'drb/drb'
gem 'rails'
#require 'rails/all'
require 'rails'
require 'google-analytics-rails'
require 'active_record'
require 'active_job'

module SushiFabric
  class Application < Rails::Application
    config.load_defaults 7.0 # for Rails 7
    # After applying Rails default settings, forcefully define action_controller
    unless config.respond_to?(:action_controller)
      config.define_singleton_method(:action_controller) do
        @action_controller ||= ActiveSupport::OrderedOptions.new
      end
    end

    # Forcefully set `config.action_controller` if it is undefined
    config.action_controller ||= ActiveSupport::OrderedOptions.new

    # default parameters
    default_root = Rails.root||Dir.pwd
    #config.workflow_manager = 'druby://localhost:12345'
    config.gstore_dir = File.join(default_root, 'public/gstore/projects')
    config.sushi_app_dir = default_root
    config.scratch_dir = '/tmp/scratch'
    config.module_source = nil
    config.course_mode = nil
    config.rails_host = nil
    config.submit_job_script_dir = nil
    config.sushi_server_class = nil
  end

  # load custmized parameters if there is
  mode = ENV['RAILS_ENV']||'production'
  config_file = File.join(Rails.root, 'config/environments', mode + ".rb")
  if File.exist?(config_file)
    SushiFabric::Application.config.action_controller ||= ActiveSupport::OrderedOptions.new
    load config_file
  else
    FileUtils.mkdir_p File.dirname(config_file)
    open(config_file+'.rb', "w") do |out|
      default_root = Rails.root||Dir.pwd
      out.print <<-EOF
module SushiFabric
  class Application < Rails::Application
    config.load_defaults 7.0 # for Rails 7
    config.action_controller ||= ActiveSupport::OrderedOptions.new
    # default parameters
    #config.workflow_manager = 'druby://localhost:12345'
    config.gstore_dir = File.join(#{default_root}, 'public/gstore/projects')
    config.sushi_app_dir = #{default_root}
    config.scratch_dir = '/tmp/scratch'
    config.module_source = nil
    config.course_mode = nil
    config.rails_host = nil
    config.submit_job_script_dir = nil
    config.sushi_server_class = nil
  end
end
      EOF
    end
  end

  config = SushiFabric::Application.config
  #WORKFLOW_MANAGER = config.workflow_manager
  SUSHI_SERVER = config.sushi_server_class
  GSTORE_DIR = config.gstore_dir
  SUSHI_APP_DIR = config.sushi_app_dir
  SCRATCH_DIR = config.scratch_dir
  MODULE_SOURCE = config.module_source
  RAILS_HOST = config.rails_host 
  SUBMIT_JOB_SCRIPT_DIR = config.submit_job_script_dir
  unless File.exist?(GSTORE_DIR)
    FileUtils.mkdir_p GSTORE_DIR
  end

  # check if there is a sqlite3 database of Ruby on Rails
  if defined?(::Project)
    NO_ROR = false
  elsif File.exist?(File.join(SUSHI_APP_DIR, "app/models")) and
    database_yml = File.join(SUSHI_APP_DIR, "config/database.yml") and
    File.exist?(database_yml)

    NO_ROR = false
    
    database_config = YAML.load(File.read(database_yml))
    db = database_config["production"]
    ActiveRecord::Base.establish_connection(
                :adapter => db["adapter"],
                :database => db["database"],
                :username => db["username"],
                :password => db["password"]
            )
    require "#{SUSHI_APP_DIR}/app/models/project"
    require "#{SUSHI_APP_DIR}/app/models/data_set"
    require "#{SUSHI_APP_DIR}/app/models/sample"
    require "#{SUSHI_APP_DIR}/app/models/job"
    require "#{SUSHI_APP_DIR}/app/jobs/application_job"
  else
    NO_ROR = true
  end

class ::Hash
  attr_reader :defaults
  alias :set :[]=
  alias :get :[]
  def []=(k1,k2,v=nil)
    if v
      @desc ||= {}
      @desc.set([k1,k2].join('_'),v)
    else
      @defaults ||= {}
      if !@defaults[k1] and k2
        if k2.instance_of?(Array)
          @defaults.set(k1,k2.first)
        elsif k2.instance_of?(Hash) and k2.first
          @defaults.set(k1,k2.first.last)
        else
          @defaults.set(k1,k2)
        end
      end
      set(k1,k2)
    end
  end
  def default_value(k,v=nil)
    if v
      @defaults[k] = v
    else
      if @defaults
        @defaults[k]
      end
    end
  end
  def data_type(k)
    if @defaults
      @defaults[k].class
    else
      v = get(k)
      if v.instance_of?(Array)
        v.first.class
      elsif v.instance_of?(Hash)
        v.values.first.class
      else
        v.class
      end
    end
  end
  def data_types
    Hash[@defaults.map{|k,v| [k, v.class]}]
  end
  def [](k1, k2=nil)
    if k2
      if @desc
        @desc.get([k1,k2].join('_'))
      else
        nil
      end
    else
      get(k1)
    end
  end
  def deep_clone_with_desc
    new_hash = self.clone
    new_hash.instance_variable_set(:@desc, @desc.clone) if @desc
    new_hash.instance_variable_set(:@defaults, @defaults.clone) if @defaults
    new_hash
  end
end
class ::String
  def tag?(tag)
    scan(/\[(.*)\]/).flatten.join =~ /#{tag}/
  end
end

class SushiApp
  attr_reader :params
  attr_reader :job_ids
  attr_reader :next_dataset_id
  attr_reader :required_columns
  attr_reader :required_params
  attr_reader :dataset_hash
  attr_reader :analysis_category
  attr_reader :description
  attr_reader :name
  attr_reader :modules
  attr_accessor :dataset_tsv_file
  attr_accessor :parameterset_tsv_file
  attr_accessor :dataset_sushi_id
  attr_accessor :data_set
  attr_accessor :project
  attr_accessor :user
  attr_accessor :name
  attr_accessor :next_dataset_name
  attr_accessor :dataset_name
  attr_accessor :next_dataset_comment
  #attr_accessor :workflow_manager
  attr_accessor :sushi_server
  attr_accessor :current_user
  attr_accessor :logger
  attr_accessor :off_bfabric_registration
  attr_accessor :mango_run_name
  attr_accessor :input_dataset_bfabric_application_number
  attr_accessor :next_dataset_bfabric_application_number
  attr_accessor :grandchild
  attr_reader :inactivate_nodes
  attr_reader :employee
  attr_accessor :queue

  def initialize
    @gstore_dir = GSTORE_DIR
    @project = nil
    @name = nil
    @params = {}
    @params['cores'] = nil
    @params['ram'] = nil
    @params['scratch'] = nil
    @params['partition'] = ''
    @params['process_mode'] = 'SAMPLE'
    @params['samples'] = ''
    @job_ids = []
    @required_columns = []
    @module_source = MODULE_SOURCE
    @modules = []
    #@workflow_manager = workflow_manager_instance||DRbObject.new_with_uri(WORKFLOW_MANAGER)
    @last_job = true
    @grandchild = true
  end
  def set_input_dataset
    if @dataset_tsv_file
      dataset_tsv = CSV.readlines(@dataset_tsv_file, :headers=>true, :col_sep=>"\t")
      @dataset_hash = []
      @dataset = []
      dataset_tsv.each do |row|
        @dataset_hash << row.to_hash
        @dataset << row.to_hash
      end

      # save in sushi db unless it is saved in sushi db
      data_set_arr = []
      headers = []
      rows = []
      dataset_name = if @dataset_name
                       @dataset_name
                     else
                       File.basename(@dataset_tsv_file).gsub(/.tsv/, '')
                     end
      data_set_arr = {'DataSetName'=>dataset_name, 'ProjectNumber'=>@project.gsub(/p/,'')}
      csv = CSV.readlines(@dataset_tsv_file, :col_sep=>"\t")
      csv.each do |row|
        if headers.empty?
          headers = row
        else
          rows << row
        end
      end
      unless NO_ROR
        @current_user ||= nil
        if @dataset_sushi_id = DataSet.save_dataset_to_database(data_set_arr: data_set_arr.to_a.flatten, headers: headers, rows: rows, user: @current_user)
          unless @off_bfabric_registration
            if dataset = DataSet.find_by_id(@dataset_sushi_id.to_i)
              dataset.register_bfabric(bfabric_application_number: @input_dataset_bfabric_application_number, update_completed_samples: true)
            end
          end
          return @dataset_sushi_id
        elsif data_set = headers[0] and data_set.instance_of?(DataSet)
          @dataset_sushi_id = data_set.id
          return @dataset_sushi_id
        end
      end
    elsif @dataset_sushi_id
      @dataset_hash = []
      @dataset = []
      if dataset = DataSet.find_by_id(@dataset_sushi_id.to_i)
        dataset.samples.each do |sample|
          @dataset_hash << sample.to_hash
          @dataset << sample.to_hash
        end
      end
      return @dataset_sushi_id
    end
    @dataset_hash
  end
  def get_columns_with_tag(tag)
    #@factor_cols = @dataset_hash.first.keys.select{|header| header =~ /\[#{tag}\]/}.map{|header| header.gsub(/\[.+\]/,'').strip}
    @dataset_hash.map{|row| 
      Hash[*row.select{|k,v| k=~/\[#{tag}\]/}.map{|k,v| [k.gsub(/\[.+\]/,'').strip,v]}.flatten]
    }
  end
  def set_default_parameters
    # this should be overwritten in a subclass
  end
  def dataset_has_column?(colname)
    flag = false
    if @dataset_hash
      @dataset_hash.map{|sample| 
        sample.each do |key, value|
          if key =~ /#{colname}/
            flag = true
          end
        end
        break
      }
    end
    flag
  end

  def set_output_files
    if @params['process_mode'] == 'SAMPLE'
      @dataset = {}
    end
    next_dataset.keys.select{|header| header.tag?('File')}.each do |header|
      @output_files ||= []
      @output_files << header
    end
    if @output_files
      @output_files = @output_files.uniq
    end
  end
  def check_required_columns
    return false unless @dataset_hash && @required_columns

    present_columns = @dataset_hash.map { |row| row.keys }.flatten.uniq
    normalized_present = present_columns.map { |col| col.gsub(/\[.+\]/, '').strip }

    # xor-mode if required_columns is nested array
    if @required_columns.all? { |entry| entry.is_a?(Array) }
      satisfied_options = @required_columns.count do |option|
        Array(option).all? { |req| normalized_present.include?(req.gsub(/\[.+\]/, '').strip) }
      end
      satisfied_options == 1
    else
      # normal and-mode
      missing = @required_columns - normalized_present
      missing.empty?
    end
  end
  def check_application_parameters
    if @required_params and (@required_params - @params.keys).empty?
      @output_params = @params.deep_clone_with_desc
      @output_params["sushi_app"] = self.class.name
      @output_params
    end
  end
  def set_user_parameters
    # this should be done in an instance of applicaiton subclass
    if @parameterset_tsv_file
      parameterset_tsv = CSV.readlines(@parameterset_tsv_file, :col_sep=>"\t")
      headers = []
      parameterset_tsv.each do |row|
        header, value = row
        headers << header
        @params[header] = if @params.data_type(header) == String or value == nil
                            value
                          else
                            eval(value)
                          end
      end
      (@params.keys - headers).each do |key|
        unless @params[key]
          @params[key] = @params.default_value(key)
        end
      end
    end
    @params
  end
  def load_project_defaults
    # Load project-specific default parameters from project_default_parametersets.tsv
    # Parameters are scoped by app name (e.g., FastqcApp::cores)
    # IMPORTANT: For Array/Hash type params (rendered as SELECT), we preserve the
    # original type and only set the 'selected' attribute to maintain UI components
    return unless @project and @gstore_dir
    
    project_defaults_file = File.join(@gstore_dir, @project, 'project_default_parametersets.tsv')
    return unless File.exist?(project_defaults_file)
    
    begin
      app_name = self.class.name
      File.readlines(project_defaults_file).each do |line|
        line.chomp!
        next if line.strip.empty?
        
        scoped_param, value = line.split("\t", 2)
        next unless scoped_param && value
        
        # Parse scoped parameter (e.g., "FastqcApp::cores")
        if scoped_param =~ /^(.+)::(.+)$/
          param_app_name = $1
          param_name = $2
          
          # Only apply parameters for the current app
          if param_app_name == app_name
            # Only set if the parameter exists in @params
            if @params.keys.include?(param_name)
              current_value = @params[param_name]
              
              # For Array/Hash types (SELECT dropdowns), preserve the component type
              # and only set the 'selected' attribute instead of overwriting the value
              if current_value.is_a?(Array) || current_value.is_a?(Hash)
                # Set the selected value without changing the component type
                @params[param_name, 'selected'] = value
              elsif current_value.is_a?(TrueClass) || current_value.is_a?(FalseClass)
                # Boolean: convert string to boolean
                @params[param_name] = (value == 'true' || value == true)
              elsif @params.data_type(param_name) == String || value.nil?
                @params[param_name] = value
              else
                # Try to eval for other types (Integer, etc.)
                begin
                  @params[param_name] = eval(value)
                rescue
                  @params[param_name] = value
                end
              end
            end
          end
        end
      end
    rescue => e
      # Silently fail if there's an error reading the file
      warn "Warning: Could not load project defaults from #{project_defaults_file}: #{e.message}"
    end
  end
  def save_project_defaults
    # Save current parameters as project-specific defaults to a temporary file
    # Returns the path to the temporary file that needs to be copied to the project folder
    # Parameters are scoped by app name (e.g., FastqcApp::cores)
    return nil unless @project and @gstore_dir
    
    project_defaults_file = File.join(@gstore_dir, @project, 'project_default_parametersets.tsv')
    
    begin
      # Read existing defaults to preserve other apps' settings
      existing_defaults = {}
      if File.exist?(project_defaults_file)
        File.readlines(project_defaults_file).each do |line|
          line.chomp!
          next if line.strip.empty?
          
          scoped_param, value = line.split("\t", 2)
          existing_defaults[scoped_param] = value if scoped_param
        end
      end
      
      # Update with current app's parameters
      app_name = self.class.name
      @params.each do |key, value|
        # Skip internal/system parameters
        next if ['process_mode', 'samples', 'node', 'sushi_app'].include?(key)
        
        scoped_param = "#{app_name}::#{key}"
        # Convert value to string representation
        value_str = if value.is_a?(String)
                      value
                    else
                      value.inspect
                    end
        existing_defaults[scoped_param] = value_str
      end
      
      # Create temporary directory with unique name to avoid conflicts
      # Using Process.pid, millisecond timestamp, and random hex for uniqueness
      # Use shared directory instead of /tmp to avoid PrivateTmp issues with Apache/Passenger
      require 'securerandom'
      shared_tmp_base = File.join(SushiFabric::Application.config.scratch_dir, "tmp")
      FileUtils.mkdir_p(shared_tmp_base) unless Dir.exist?(shared_tmp_base)
      temp_dir = File.join(shared_tmp_base, "project_defaults_#{Process.pid}_#{Time.now.strftime("%Y%m%d%H%M%S%L")}_#{SecureRandom.hex(4)}")
      FileUtils.mkdir_p(temp_dir)
      
      # Write to file with the correct final name
      temp_file = File.join(temp_dir, 'project_default_parametersets.tsv')
      File.open(temp_file, 'w') do |f|
        existing_defaults.sort.each do |scoped_param, value|
          f.puts "#{scoped_param}\t#{value}"
        end
      end
      
      # Set appropriate permissions
      FileUtils.chmod(0664, temp_file) if File.exist?(temp_file)
      
      return temp_file
    rescue => e
      # Log error but don't fail the job
      warn "Warning: Could not save project defaults to temporary file: #{e.message}"
      return nil
    end
  end
  def set_dir_paths
    ## sushi figures out where to put the resulting dataset
    unless @name and @project
      raise "should set #name and #project"
    end
    @name.gsub!(/\s/,'_')
    @result_dir_base = if @next_dataset_name
                        [@next_dataset_name, Time.now.strftime("%Y-%m-%d--%H-%M-%S")].join("_")
                      else
                        [@name, @dataset_sushi_id.to_s, Time.now.strftime("%Y-%m-%d--%H-%M-%S")].join("_")
                      end
    @result_dir = File.join(@project, @result_dir_base)
    @scratch_result_dir = File.join(SCRATCH_DIR, @result_dir_base)
    @job_script_dir = File.join(@scratch_result_dir, 'scripts')
    @uploaded_files_dir = File.join(@scratch_result_dir, 'uploaded')
    @gstore_result_dir = File.join(@gstore_dir, @result_dir)
    @gstore_script_dir = File.join(@gstore_result_dir, 'scripts')
    @gstore_project_dir = File.join(@gstore_dir, @project)
    @gstore_uploaded_dir = File.join(@gstore_result_dir, 'uploaded')
    set_file_paths
  end
  def prepare_result_dir
    FileUtils.mkdir_p(@scratch_result_dir)
    FileUtils.mkdir_p(@job_script_dir)
    @uploaded_files = []
    @params.each do |key, value|
      if @params[key, 'file_upload']
        FileUtils.mkdir_p(@uploaded_files_dir)
        @uploaded_files << value
      end
    end
  end
  def check_latest_modules_version(modules)
    command_out =  %x[ bash -lc "source #{@module_source}; module whatis #{modules.join(" ")} 2>&1" | cut -f 1 -d " " | uniq ]
    latest_modules = []
    command_out.split("\n").each do |line_|
      line = line_.chomp
      unless line.empty?
        if line =~ /#{modules.join("|")}/
          latest_modules << line
        end
      end
    end
    latest_modules
  end
  def job_header
    @scratch_dir = if @params['process_mode'] == 'SAMPLE'
                     @scratch_result_dir + "_" + @dataset['Name'] + '_temp$$'
                   else
                     @scratch_result_dir + '_temp$$'
                   end
    hold_jid_option = if @dataset_sushi_id and parent_data_set = DataSet.find_by_id(@dataset_sushi_id.to_i) and !parent_data_set.jobs.empty? and parent_data_set_job_ids = parent_data_set.jobs.map{|job| job.submit_job_id} and !parent_data_set_job_ids.join.empty?
                        "#SBATCH --dependency=afterany:#{parent_data_set_job_ids.join(":")}"
                      else
                        ''
                      end
    module_src_command = if @module_source and @modules and !@modules.empty?
                       "source #{@module_source}"
                     else
                       ""
                     end
    module_add_commands = if @modules and !@modules.empty?
                            modules_with_version = check_latest_modules_version(@modules)
                            if @modules.length == modules_with_version.length
                              modules_with_version.compact!
                              "module add #{modules_with_version.join(' ')}"
                            else
                              @logger.error("#"*100)
                              @logger.error("# Error in checking modules ")
                              @logger.error("# Please check if all modules are correctly installed, searched #{@modules.join(",")} but only detected #{modules_with_version.join(",")}")
                              @logger.error("#"*100)
                              # "exit # Please check if all modules are correctly installed, searched #{@modules.join(",")} but only detected #{modules_with_version.join(",")}"
                              ""
                            end
                          else
                            ""
                          end
    @out.print <<-EOF
#!/bin/bash
#{hold_jid_option}
set -eux
set -o pipefail
umask 0002

#### SET THE STAGE
SCRATCH_DIR=#{@scratch_dir}
GSTORE_DIR=#{@gstore_dir}
INPUT_DATASET=#{@input_dataset_tsv_path}
LAST_JOB=#{@last_job.to_s.upcase}
echo "Job runs on `hostname`"
echo "at $SCRATCH_DIR"
mkdir $SCRATCH_DIR || exit 1
cd $SCRATCH_DIR || exit 1
#{module_src_command}
#{module_add_commands}

    EOF
  end
  def job_footer
    @out.print "#### JOB IS DONE WE PUT THINGS IN PLACE AND CLEAN AUP\n"
    if File.exist?("/usr/local/ngseq/miniconda3/etc/profile.d/conda.sh")
      @out.print <<-EOS
. "/usr/local/ngseq/miniconda3/etc/profile.d/conda.sh"
conda activate sushi
      EOS
    end

    src_files = []
    dest_dirs = []
    greq = (copy_commands("AAA", "BBB").join =~ /g-req/)
    if @output_files
      @output_files.map{|header| next_dataset[header]}.each do |file|
        src_file = File.basename(file)
        dest_dir = File.dirname(File.join(@gstore_dir, file))
        src_files << src_file
        dest_dirs << dest_dir
      end
      if dest_dirs.uniq.length == 1 and greq
        src_file = src_files.join(" ")
        dest_dir = dest_dirs.first
        @out.print copy_commands(src_file, dest_dir, nil, @queue).join("\n"), "\n"
      else
        @output_files.map{|header| next_dataset[header]}.each do |file|
          # in actual case, to save under /srv/gstore/
          src_file = File.basename(file)
          dest_dir = File.dirname(File.join(@gstore_dir, file))
          @out.print copy_commands(src_file, dest_dir, nil, @queue).join("\n"), "\n"
        end
      end
    end
    # Copy grandchild dataset files
    grandchild_datasets.each do |ds|
      ds.keys.select{|header| header.tag?('File')}.each do |header|
        file = ds[header]
        next if file.to_s.empty?
        src_file = File.basename(file)
        dest_dir = File.dirname(File.join(@gstore_dir, file))
        @out.print copy_commands(src_file, dest_dir, nil, @queue).join("\n"), "\n"
      end
    end
    @out.print <<-EOF
cd #{SCRATCH_DIR}
rm -rf #{@scratch_dir} || exit 1

    EOF

  end
  def job_main
    @out.print "#### NOW THE ACTUAL JOBS STARTS\n"
    @out.print commands, "\n\n"
  end
  def next_dataset
    # this should be overwritten in a subclass
  end
  def grandchild_datasets
    # this should be overwritten in a subclass
    # returns an array of hashes
    []
  end
  def commands
    # this should be overwritten in a subclass
  end
  def generate_new_job_script(script_name, script_content)
    new_job_script = File.basename(script_name) + "_" + Time.now.strftime("%Y%m%d%H%M%S%L")
    new_job_script = File.join(SUBMIT_JOB_SCRIPT_DIR, new_job_script)
    open(new_job_script, 'w') do |out|
      out.print script_content
      out.print "\necho __SCRIPT END__\n"
    end
    new_job_script
  end
  def submit_job_command(script_file, script_content, option='')
    if script_name = File.basename(script_file) and script_name =~ /\.sh/
      script_name = script_name.split(/\.sh/).first + ".sh"
      new_job_script = generate_new_job_script(script_name, script_content)
      new_job_script_base = File.basename(new_job_script)
      log_file = File.join(SUBMIT_JOB_SCRIPT_DIR, new_job_script_base + "_o.log")
      err_file = File.join(SUBMIT_JOB_SCRIPT_DIR, new_job_script_base + "_e.log")
      command = "sbatch -o #{log_file} -e #{err_file} -N 1 #{option} #{new_job_script}"
      [command, new_job_script, log_file, err_file]
    else
      err_msg = "submit_job_command, ERROR: script_name is not *.sh: #{File.basename(script_file)}"
      warn err_msg
      raise err_msg
    end
  end
  def submit(script_path, mock=false)
    sbatch_options = []
    sbatch_options << "--mem=#{@params['ram']}G" unless @params['ram'].to_s.empty?
    sbatch_options << "-n #{@params['cores']}" unless @params['cores'].to_s.empty?
    sbatch_options << "--gres=scratch:#{@params['scratch']}" unless @params['scratch'].to_s.empty?
    sbatch_options << "-p #{@params['partition']}" unless @params['partition'].to_s.empty?
    sbatch_options << "--nice=#{@params['nice']}" unless @params['nice'].to_s.empty?

    script_content = File.read(script_path)

    submit_command, new_script_path, stdout_path, stderr_path = submit_job_command(script_path, script_content, sbatch_options.join(' '))
    [submit_command, new_script_path, stdout_path, stderr_path]
  end
  def preprocess
    # this should be overwritten in a subclass
  end
  def set_file_paths
    @parameter_file = 'parameters.tsv'
    @input_dataset_file = 'input_dataset.tsv'
    @next_dataset_file = 'dataset.tsv'
    @input_dataset_tsv_path = File.join(@gstore_result_dir, @input_dataset_file)
    @parameters_tsv_path = File.join(@gstore_result_dir, @parameter_file)
    @next_dataset_tsv_path = File.join(@gstore_result_dir, @next_dataset_file)
  end
  def save_parameters_as_tsv
    file_path = File.join(@scratch_result_dir, @parameter_file)
    CSV.open(file_path, 'w', :col_sep=>"\t") do |out|
      @output_params.each do |key, value|
        if @output_params[key, 'file_upload'] and !value.to_s.empty?
          uploaded_file_path = File.join(@result_dir, "uploaded", File.basename(value))
          out << [key, uploaded_file_path]
          @params[key] = uploaded_file_path
          @output_params[key] = uploaded_file_path
        else
          out << [key, value]
        end
      end
    end
    file_path
  end
  def save_input_dataset_as_tsv
    file_path = File.join(@scratch_result_dir, @input_dataset_file)
    CSV.open(file_path, 'w', :col_sep=>"\t") do |out|
      headers = @dataset_hash.map{|row| row.keys}.flatten.uniq
      out << headers
      @dataset_hash.each do |row|
        out << headers.map{|header| 
          val = row[header]
          val.to_s.empty? ? nil:val
        }
      end
    end
    file_path
  end
  def save_next_dataset_as_tsv
    headers = @result_dataset.map{|row| row.keys}.flatten.uniq
    file_path = File.join(@scratch_result_dir, @next_dataset_file)
    CSV.open(file_path, 'w', :col_sep=>"\t") do |out|
      out << headers
      @result_dataset.each do |row_hash|
        out << headers.map{|header| 
          val = row_hash[header]
          val.to_s.empty? ? nil:val
        }
      end
    end
    file_path
  end
  def save_grandchild_datasets_as_tsv
    grandchild_data = grandchild_datasets
    return [] if grandchild_data.empty?

    grandchild_tsv_paths = []
    ## save multi-sample dataset
    headers = grandchild_data.map{|row| row.keys}.flatten.uniq
    file_path = File.join(@scratch_result_dir, "grandchild_dataset.tsv")
    CSV.open(file_path, 'w', :col_sep=>"\t") do |out|
      out << headers
      grandchild_data.each do |row_hash|
        out << headers.map{|header| 
          val = row_hash[header]
          val.to_s.empty? ? nil:val
        }
      end

#    grandchild_data.each_with_index do |dataset_hash, index|
#      file_name = "grandchild_dataset_#{index + 1}.tsv"
#      file_path = File.join(@scratch_result_dir, file_name)
      
#      CSV.open(file_path, 'w', :col_sep=>"\t") do |out|
        # Write header
#        out << dataset_hash.keys
        # Write single row (since each grandchild_dataset represents one dataset)
#        out << dataset_hash.keys.map { |key| 
#          val = dataset_hash[key]
#          val.to_s.empty? ? nil : val
#        }
#      end
      
      grandchild_tsv_paths << file_path
    end
    
    grandchild_tsv_paths
  end
  def copy_commands(org_dir, dest_parent_dir, now=nil, queue="light")
    @sushi_server||=eval(SushiFabric::Application.config.sushi_server_class).new
    com = ''
    cnt_retry = 0
    begin
      com = @sushi_server.copy_commands(org_dir, dest_parent_dir, now, queue)
    rescue => e
      time = Time.now.strftime("[%Y.%m.%d %H:%M:%S]")
      @logger.error("*"*50)
      @logger.error("copy_command error #{time}")
      @logger.error("error: #{e}")
      @logger.error("org_dir: #{org_dir}, dest_parent_dir: #{dest_parent_dir}, now: #{now}")
      @logger.error("*"*50)
      sleep 1
      cnt_retry += 1
      retry if cnt_retry < 3
    end
    com
  end
  def copy_uploaded_files
    if not @uploaded_files.empty?
      @uploaded_files.compact.select{|file| !file.empty?}.each do |file|
        FileUtils.cp(file, @uploaded_files_dir)
        command = "cp #{file} #{@uploaded_files_dir}"
        warn command
        FileUtils.rm_r(File.dirname(file))
        command = "rm -rf #{File.dirname(file)}"
        warn command
      end
    end
  end
  def copy_inputdataset_parameter_jobscripts
    org = @scratch_result_dir
    dest = @gstore_project_dir
    copy_commands(org, dest, 'now').each do |command|
      warn `which python`
      warn command
      unless system command
        raise "fails in copying input_dataset, parameters and jobscript files from /scratch to /gstore"
      end
    end
    #sleep 1
  end
  def copy_nextdataset
    org = @next_dataset_tsv_path
    dest = File.join(@gstore_project_dir, @result_dir_base)
    copy_commands(org, dest, 'now').each do |command|
      warn `which python`
      warn command
      unless system command
        raise "fails in copying next_dataset files from /scratch to /gstore"
      end
    end
    sleep 1
    command = "rm -rf #{@scratch_result_dir}"
    `#{command}`
  end
  #def cluster_nodes
  #  @workflow_manager||=DRbObject.new_with_uri(WORKFLOW_MANAGER)
  #  @workflow_manager.cluster_nodes
  #end
  #def default_node
  #  @workflow_manager||=DRbObject.new_with_uri(WORKFLOW_MANAGER)
  #  @workflow_manager.default_node
  #end

  def make_job_script(append = false)
    @out = if append
             open(@job_script, 'a')
           else
             open(@job_script, 'w')
           end
    job_header
    job_main
    job_footer
    @out.close
  end
  def sample_mode
    selected_samples = unless @params['samples'].empty?
                         Hash[*@params['samples'].split(',').map{|sample_name| [sample_name, true]}.flatten]
                       else
                         Hash[*@dataset_hash.map{|row| row['Name']}.map{|sample_name| [sample_name, true]}.flatten]
                       end
    @dataset_hash.each_with_index do |row, i|
      @dataset = Hash[*row.map{|key,value| [key.gsub(/\[.+\]/,'').strip, value]}.flatten]
      if selected_samples[@dataset['Name']]
        sample_name = @dataset['Name']||@dataset.first
        @job_script = if @dataset_sushi_id and dataset = DataSet.find_by_id(@dataset_sushi_id.to_i)
                        File.join(@job_script_dir, @analysis_category + '_' + sample_name) + '_' + dataset.name.gsub(/\s+/,'_') + "_" + @dataset_sushi_id.to_s + '.sh'
                      else
                        File.join(@job_script_dir, @analysis_category + '_' + sample_name) + '.sh'
                      end
        @last_job = (i == @dataset_hash.length - 1)
        make_job_script
        @job_scripts << @job_script
        @result_dataset << next_dataset
      end
    end
  end
  def dataset_mode
    selected_samples = unless @params['samples'].empty?
                         Hash[*@params['samples'].split(',').map{|sample_name| [sample_name, true]}.flatten]
                       else
                         Hash[*@dataset_hash.map{|row| row['Name']}.map{|sample_name| [sample_name, true]}.flatten]
                       end
    # for a case of @dataset is used in def next_datast in SUSHIApp
    @dataset = []
    @dataset_hash.each do |row|
      if selected_samples[row['Name']]
        @dataset << row
      end
    end
    @job_script = if @dataset_sushi_id and dataset = DataSet.find_by_id(@dataset_sushi_id.to_i)
                    File.join(@job_script_dir, @analysis_category + '_' + dataset.name.gsub(/[\s+,\/]/,'_') + "_" + @dataset_sushi_id.to_s + '.sh')
                  else 
                    File.join(@job_script_dir, @analysis_category + '_' + 'job_script.sh')
                  end
    make_job_script
    @job_scripts << @job_script
    @result_dataset << next_dataset
  end
  def batch_mode
    @job_script = if @dataset_sushi_id and dataset = DataSet.find_by_id(@dataset_sushi_id.to_i)
                    File.join(@job_script_dir, dataset.name.gsub(/\s+/,'_') + '.sh')
                  else 
                    File.join(@job_script_dir, 'job_script.sh')
                  end
    @dataset_hash.each do |row|
      @dataset = Hash[*row.map{|key,value| [key.gsub(/\[.+\]/,'').strip, value]}.flatten]
      make_job_script('append')
      @result_dataset << next_dataset
    end
    @job_scripts << @job_script
  end
  def save_parameters_in_sushi_db
    if @next_dataset_id and next_dataset = DataSet.find_by_id(@next_dataset_id)
      next_dataset.job_parameters = @output_params
      next_dataset.save
    end
  end
  def save_grandchild_datasets_to_database
    return if @grandchild_dataset_tsv_paths.nil? || @grandchild_dataset_tsv_paths.empty? || !@next_dataset_id

    grandchild_data = grandchild_datasets
    @grandchild_dataset_ids = []

    @grandchild_dataset_tsv_paths.each_with_index do |tsv_path, index|
      headers = []
      rows = []
      
      csv = CSV.readlines(tsv_path, :col_sep=>"\t")
      csv.each do |row|
        if headers.empty?
          headers = row
        else
          rows << row
        end
      end
      
      grandchild_dataset_name = if @params && @params['grandchildName'] && !@params['grandchildName'].to_s.strip.empty?
                                  @params['grandchildName']
                                elsif grandchild_data[index] && grandchild_data[index]['Name'] && !grandchild_data[index]['Name'].to_s.strip.empty?
                                  grandchild_data[index]['Name']
                                else
                                  "#{@name}_grandchild_#{index + 1}"
                                end
      grandchild_comment = "autogenerated grandchild" #"Grandchild dataset #{index + 1}"
      
      data_set_arr = {
        'DataSetName' => grandchild_dataset_name, 
        'ProjectNumber' => @project.gsub(/p/,''), 
        'ParentID' => @next_dataset_id,
        'Comment' => grandchild_comment.to_s
      }
      
      unless NO_ROR
        # Determine child flag: if @grandchild is true, force true; else inherit @child
        grandchild_child_flag = @grandchild ? true : !!@child
        grandchild_dataset_id = DataSet.save_dataset_to_database(
          data_set_arr: data_set_arr.to_a.flatten, 
          headers: headers, 
          rows: rows, 
          user: @current_user, 
          child: grandchild_child_flag, 
          sushi_app_name: self.class.name
        )
        @grandchild_dataset_ids << grandchild_dataset_id
      end
    end
  end
  def main(mock=false)
    ## sushi writes creates the job scripts and builds the result data set that is to be generated
    @result_dataset = []
    @job_scripts = []
    if @params['process_mode'] == 'SAMPLE'
      sample_mode
    elsif @params['process_mode'] == 'DATASET'
      dataset_mode
    elsif @params['process_mode'] == 'BATCH'
      batch_mode
    else 
      #stop
      warn "the process mode (#{@params['process_mode']}) is not defined"
      raise "stop job submitting"
    end
    if mock
      make_dummy_files
    end

    warn 'result dataset: '
    warn @result_dataset

    # copy application data to gstore 
    @next_dataset_tsv_path = save_next_dataset_as_tsv
    
    # Save grandchild datasets as TSV files
    @grandchild_dataset_tsv_paths = save_grandchild_datasets_as_tsv

    if @dataset_sushi_id and dataset = DataSet.find_by_id(@dataset_sushi_id.to_i)
      @project_id = dataset.project.id
      data_set_arr = []
      headers = []
      rows = []
      next_dataset_name = if name = @next_dataset_name
                            name.to_s
                          else
                            "#{@name.gsub(/\s/,'').gsub(/_/,'')}_#{dataset.id}"
                          end
      data_set_arr = if @params['next_dataset_root']
                       {'DataSetName'=>next_dataset_name, 'ProjectNumber'=>@project.gsub(/p/,''), 'Comment'=>@next_dataset_comment.to_s}
                     else
                       {'DataSetName'=>next_dataset_name, 'ProjectNumber'=>@project.gsub(/p/,''), 'ParentID'=>@dataset_sushi_id, 'Comment'=>@next_dataset_comment.to_s}
                     end

      csv = CSV.readlines(@next_dataset_tsv_path, :col_sep=>"\t")
      csv.each do |row|
        if headers.empty?
          headers = row
        else
          rows << row
        end
      end
      unless NO_ROR
        @current_user ||= nil
        @next_dataset_id = DataSet.save_dataset_to_database(data_set_arr: data_set_arr.to_a.flatten, headers: headers, rows: rows, user: @current_user, child: @child, sushi_app_name: self.class.name)
        save_parameters_in_sushi_db
        
        # Save grandchild datasets to database
        save_grandchild_datasets_to_database
      end
    end
    copy_uploaded_files
    copy_inputdataset_parameter_jobscripts

    # job submittion
    gstore_job_script_paths = []
    #submit_jobs = []
    @job_scripts.each_with_index do |job_script, i|
      #submit_command, script_path, stdout_path, stderr_path = submit(job_script, mock)
      #submit_jobs << [submit_command, script_path, stdout_path, stderr_path]
      gstore_job_script_paths << File.join(@gstore_script_dir, File.basename(job_script))
    end

    warn ""
    warn 'job scripts: '
    warn @job_scripts
    warn 'gstore_job_script_paths: '
    warn gstore_job_script_paths

    #unless @job_ids.empty? or NO_ROR
    #unless submit_jobs.empty?
    unless gstore_job_script_paths.empty?
      # save job and dataset relation in Sushi DB
      gstore_job_script_paths.each do |script_path|
        unless mock
          new_job = Job.new
          new_job.script_path = script_path
          new_job.next_dataset_id = @next_dataset_id
          new_job.status = "CREATED"
          new_job.user = (@user || "sushi_lover")
          new_job.input_dataset_id = dataset.id if dataset
          new_job.save
          #if dataset
          #  new_job.data_set = dataset
          #  new_job.data_set.jobs << new_job
          #  new_job.data_set.save
          #end
        end
      end
    end

    copy_nextdataset
  end
  def run
    dataset_id = test_run

    ## the user presses RUN
    prepare_result_dir

    ## copy application data to gstore 
    save_parameters_as_tsv
    save_input_dataset_as_tsv

    main
    dataset_id
  end
  def make_dummy_files
    dummy_files_header = []
    headers = @result_dataset.map{|row| row.keys}.flatten.uniq
    headers.select{|header| header.tag?('File')||header.tag?('Link')}.each do |header|
      dummy_files_header << header
    end
    dummy_files_ = []
    @result_dataset.each do |row|
      dummy_files_.concat(dummy_files_header.map{|header| row[header]})
    end
    dummy_files = []
    dummy_files_.each do |file|
      dummy_files << file.gsub(@result_dir, '')
    end
    dummy_files.uniq!

    dirs = []
    dummy_files.permutation(2).each do |a,b|
      if a.include?(b) and b !~ /\./
        dirs << b
      end
    end
    dirs.each do |dir|
      dummy_files.delete(dir)
    end
    dirs.each do |dir|
      command = "mkdir -p #{File.join(@scratch_result_dir, dir)}"
      warn command
      `#{command}`
    end
    dummy_files.each do |file|
      command = if file =~ /.html/
                  "echo 'Hello, SUSHI world!' > #{File.join(@scratch_result_dir, file)}"
                else
                  "touch #{File.join(@scratch_result_dir, file)}"
                end
      warn command
      `#{command}`
    end
  end
  def mock_run
    test_run
    prepare_result_dir
    save_parameters_as_tsv
    save_input_dataset_as_tsv
    main(true)
  end
  def test_run
    dataset_id = set_input_dataset
    set_dir_paths
    preprocess
    set_output_files
    set_user_parameters

    failures = 0
    err_msgs = []
    warn 'check project name: '
    unless @project
      err_msg = []
      err_msg << "\e[31mFAILURE\e[0m: project number is required but not found. you should set it in usecase."
      err_msg << "\tex.)"
      err_msg << "\tapp = #{self.class}.new"
      err_msg << "\tapp.project = 'p1001'"
      warn err_msg.join("\n")
      err_msgs.concat(err_msg)
      failures += 1
    else
      warn "\e[32mPASSED\e[0m:\n\t@project=#{@project}"
    end

    warn 'check user name: '
    unless @user
      err_msg = []
      err_msg << "\e[31mWARNING\e[0m: user number is ought to be added but not found. you should set it in usecase. Default will be 'sushi lover'"
      err_msg << "\tex.)"
      err_msg << "\tapp = #{self.class}.new"
      err_msg << "\tapp.user = 'masa'"
      warn err_msg.join("\n")
      err_msgs.concat(err_msg)
    else
      warn "\e[32mPASSED\e[0m:\n\t@user=#{@user}"
    end

    warn 'check application name: '
    if @name.to_s.empty?
      err_msg = []
      err_msg << "\e[31mFAILURE\e[0m: application name is required but not found. you should set it in application class."
      err_msg << "\tex.)"
      err_msg << "\tclass #{self.class}"
      err_msg << "\t def initialize"
      err_msg << "\t  @name = '#{self.class}'"
      err_msg << "\t end"
      err_msg << "\tend"
      warn err_msg.join("\n")
      err_msgs.concat(err_msg)
      failures += 1
    else
      warn "\e[32mPASSED\e[0m:\n\t@name=#{@name}"
    end

    warn 'check analysis_category: '
    if @analysis_category.to_s.empty?
      err_msg = []
      err_msg << "\e[31mFAILURE\e[0m: analysis_category is required but not found. you should set it in application class."
      err_msg << "\tex.)"
      err_msg << "\tclass #{self.class}"
      err_msg << "\t def initialize"
      err_msg << "\t  @analysis_category = 'Mapping'"
      err_msg << "\t end"
      err_msg << "\tend"
      warn err_msg.join("\n")
      err_msgs.concat(err_msg)
      failures += 1
    else
      warn "\e[32mPASSED\e[0m:\n\t@analysis_category=#{@analysis_category}"
    end

    warn 'check dataset: '
    if !@dataset_hash or @dataset_hash.empty?
      err_msg = []
      err_msg << "\e[31mFAILURE\e[0m: dataset is not found. you should set it by using #{self.class}#dataset_sushi_id or #{self.class}#dataset_tsv_file properties"
      err_msg << "\tex.)"
      err_msg << "\tusecase = #{self.class}.new"
      err_msg << "\tusecase.dataset_tsv_file = \"dataset.tsv\""
      warn err_msg.join("\n")
      err_msgs.concat(err_msg)
      failures += 1
    else
      warn "\e[32mPASSED\e[0m:\n\t@dataset_hash.length = #{@dataset_hash.length}"
    end

    warn 'check required columns: '
    unless check_required_columns
      err_msg = []
      err_msg << "\e[31mFAILURE\e[0m: required_column(s) is not found in dataset. you should set it in application class."
      err_msg << "\tex.)"
      err_msg << "\tdef initialize"
      err_msg << "\t  @required_columns = ['Name', 'Read1']"
      warn err_msg.join("\n")
      err_msgs.concat(err_msg)
      failures += 1
    else
      warn "\e[32mPASSED\e[0m:"
    end
    warn "\trequired columns: #{@required_columns}"
    warn "\tdataset  columns: #{@dataset_hash.map{|row| row.keys}.flatten.uniq}" if @dataset_hash

    warn 'check required parameters: '
    unless check_application_parameters
      err_msg = []
      err_msg << "\e[31mFAILURE\e[0m: required_param(s) is not set yet. you should set it in usecase"
      err_msg << "\tmissing params: #{@required_params-@params.keys}" if @required_params
      err_msg << "\tex.)"
      err_msg << "\tusecase = #{self.class}.new"
      if @required_params
        err_msg << "\tusecase.params['#{(@required_params-@params.keys)[0]}'] = parameter"
      else
        err_msg << "\tusecase.params['parameter name'] = default_parameter"
      end
      warn err_msg.join("\n")
      err_msgs.concat(err_msg)
      failures += 1
    else
      warn "\e[32mPASSED\e[0m:"
    end
    warn "\tparameters: #{@params.keys}"
    warn "\trequired  : #{@required_params}"

    warn 'check next dataset: '
    if @params['process_mode'] == 'SAMPLE'
      @dataset={}
    end
    unless self.next_dataset
      err_msg = []
      err_msg << "\e[31mFAILURE\e[0m: next dataset is not set yet. you should overwrite SushiApp#next_dataset method in #{self.class}"
      err_msg << "\tnote: the return value should be Hash (key: column title, value: value in a tsv table)"
      warn err_msg.join("\n")
      err_msgs.concat(err_msg)
      failures += 1
    else
      warn "\e[32mPASSED\e[0m:"
    end

    # Add grandchild datasets check
    warn 'check grandchild datasets: '
    grandchild_data = self.grandchild_datasets
    if grandchild_data.any?
      # Validate that each grandchild dataset has a unique Name
      names = grandchild_data.map { |ds| ds['Name'] }.compact
      if names.length != names.uniq.length
        err_msg = []
        err_msg << "\e[31mFAILURE\e[0m: grandchild datasets must have unique 'Name' values"
        err_msg << "\tfound names: #{names}"
        warn err_msg.join("\n")
        err_msgs.concat(err_msg)
        failures += 1
      else
        warn "\e[32mPASSED\e[0m: #{grandchild_data.length} grandchild dataset(s) defined"
      end
    else
      warn "\e[32mINFO\e[0m: No grandchild datasets defined"
    end

    warn 'check output files: '
    if !@output_files or @output_files.empty?
      err_msg = []
      err_msg << "\e[31mWARNING\e[0m: no output files. you will not get any output files after the job running. you can set @output_files (array) in #{self.class}"
      err_msg << "\tnote: usually it should be define in initialize method"
      err_msg << "\t      the elements of @output_files should be chosen from #{self.class}#next_dataset.keys"
      err_msg << "\t      #{self.class}#next_dataset.keys: #{self.next_dataset.keys}" if self.next_dataset
      warn err_msg.join("\n")
      err_msgs.concat(err_msg)
    else
      warn "\e[32mPASSED\e[0m:"
    end

    warn 'check commands: '
    if @params['process_mode'] == 'SAMPLE'
      @dataset_hash.each do |row|
        @dataset = Hash[*row.map{|key,value| [key.gsub(/\[.+\]/,'').strip, value]}.flatten]
        unless com = commands
          err_msg = []
          err_msg << "\e[31mFAILURE\e[0m: any commands is not defined yet. you should overwrite SushiApp#commands method in #{self.class}"
          err_msg << "\tnote: the return value should be String (this will be in the main body of submitted job script)"
          warn err_msg.join("\n")
          err_msgs.concat(err_msg)
          failures += 1
        else
          warn "\e[32mPASSED\e[0m:"
          warn "generated command will be:"
          warn "\t"+com.split(/\n/).join("\n\t")+"\n"
        end
      end
    elsif @params['process_mode'] == 'DATASET'
      unless com = commands
        err_msg = []
        err_msg << "\e[31mFAILURE\e[0m: any commands is not defined yet. you should overwrite SushiApp#commands method in #{self.class}"
        err_msg << "\tnote: the return value should be String (this will be in the main body of submitted job script)"
        warn err_msg.join("\n")
        err_msgs.concat(err_msg)
        failures += 1
      else
        warn "\e[32mPASSED\e[0m:"
        warn "generated command will be:"
        warn "\t"+com.split(/\n/).join("\n\t")+"\n"
      end
    end

    warn 'check sushi server instance: '
    begin
      #@workflow_manager||=DRbObject.new_with_uri(WORKFLOW_MANAGER)
      @sushi_server||=eval(SushiFabric::Application.config.sushi_server_class).new
      hello = @sushi_server.hello
    rescue
    end
    unless hello =~ /hello/
      err_msg = "\e[31mFAILURE\e[0m: sushi server does not reply. check if sushi server is working"
      warn err_msg
      err_msgs.concat([err_msg])
      failures += 1
    else
      #warn "\e[32mPASSED\e[0m: #{WORKFLOW_MANAGER}"
      warn "\e[32mPASSED\e[0m: #{SUSHI_SERVER}"
    end

    if failures > 0
      warn ""
      err_msg = "\e[31mFailures (#{failures})\e[0m: All failures should be solved"
      warn err_msg
      err_msgs.unshift(err_msg)
      raise "\n"+err_msgs.join("\n")+"\n\n"
    else
      warn "All checks \e[32mPASSED\e[0m"
    end

    dataset_id
  end
end

class SushiServerClass
    def copy_commands(org_dir, dest_parent_dir, now=nil, queue="light")
    end
    def delete_command(target)
    end
    def hello
      "hello"
    end
    def command_available?(cmd)
      system("which #{cmd} > /dev/null 2>&1")
    end
end

class DemoSushi < SushiServerClass
  def copy_commands(org_dir, dest_parent_dir, now=nil, queue="light")
    commands = ["rsync -r #{org_dir} #{dest_parent_dir}/"]
  end
  def delete_command(target, bfabric_dataset_id: nil)
    command = "rm -rf #{target}"
  end
end

CourseSushi = DemoSushi

class ProdSushi < SushiServerClass
  def copy_commands(org_dir, dest_parent_dir, now=nil, queue="light")
    commands = if now == "force"
                 target_file = File.join(dest_parent_dir, File.basename(org_dir))
                 ["g-req copynow -f #{org_dir} #{dest_parent_dir}"]
               elsif now
                 ["g-req copynow #{org_dir} #{dest_parent_dir}"]
               elsif queue.nil? or queue == "light"
                 ["g-req -w copy #{org_dir} #{dest_parent_dir}"]
               else
                 ["g-req -w copy -f heavy #{org_dir} #{dest_parent_dir}"]
               end
  end
  def delete_command(target, bfabric_dataset_id: nil)
    command = if command_available?("delete_workunit_by_dataset_id") and bfabric_dataset_id
                "g-req remove --bfabric-dataset-id #{bfabric_dataset_id} #{target}"
              else
                "g-req remove #{target}"
              end
  end
end

TestSushi = ProdSushi

end # SushiFabric
