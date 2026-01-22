require 'yaml'
require 'fileutils'
require_relative 'nf_core_info_fetcher'
require_relative 'global_variables'

module NfCoreAppFactory
  PIPELINE_CONFIG_PATH = File.expand_path('../../config/nf_core_pipelines.yml', __FILE__)
  LIB_DIR = File.expand_path('..', __FILE__)
  
  # Minimal required_columns - only 'Name' is required for maximum compatibility
  # nf-core pipelines have varied input requirements, and users should ensure
  # their dataset has appropriate columns for the selected pipeline.
  # Specific requirements can be overridden in nf_core_pipelines.yml
  DEFAULTS = {
    'analysis_category' => 'Misc',
    'required_columns' => ['Name'],
    'params' => { 'cores' => 8, 'ram' => 30, 'scratch' => 100 }
  }
  
  CATEGORY_MAPPING = {
    'rna-seq' => 'Transcriptomics',
    'dna-seq' => 'Genomics',
    'chip-seq' => 'Epigenetics',
    'metagenomics' => 'Metagenomics',
    'quality-control' => 'QC',
    'sc-rna-seq' => 'SingleCell',
    'single-cell' => 'SingleCell'
  }
  
  def self.load_pipelines
    # Get all pipelines with basic info from single API call (efficient)
    all_pipelines = NfCoreInfoFetcher.fetch_all_pipelines_with_info || {}
    puts "NfCoreAppFactory: API returned #{all_pipelines.size} pipelines"
    $stdout.flush
    
    # Load manual config and defaults
    manual_config = {}
    defaults_config = {}
    if File.exist?(PIPELINE_CONFIG_PATH)
      yaml_content = YAML.load_file(PIPELINE_CONFIG_PATH) || {}
      manual_config = yaml_content['pipelines'] || {}
      defaults_config = yaml_content['defaults'] || {}
      puts "NfCoreAppFactory: Manual config has #{manual_config.size} pipelines"
      $stdout.flush
    end
    
    # Deep merge manual config to override API values
    manual_config.each do |key, config|
      if all_pipelines[key]
        all_pipelines[key] = deep_merge(all_pipelines[key], config)
      else
        all_pipelines[key] = config.merge({'nf_core_name' => key})
      end
    end
    
    # Apply defaults to all pipelines that don't have explicit config
    all_pipelines.each do |key, pipeline_info|
      # Apply default params if not set
      unless pipeline_info['params']
        pipeline_info['params'] = defaults_config['params'] || DEFAULTS['params']
      end
      
      # Apply default samplesheet_mapping if not set and input_type is samplesheet
      unless pipeline_info['samplesheet_mapping']
        pipeline_info['samplesheet_mapping'] = defaults_config['samplesheet_mapping']
      end
    end
    
    puts "NfCoreAppFactory: Total pipelines to register: #{all_pipelines.size}"
    $stdout.flush
    all_pipelines
  end
  
  # Load defaults from config file
  def self.load_defaults
    if File.exist?(PIPELINE_CONFIG_PATH)
      yaml_content = YAML.load_file(PIPELINE_CONFIG_PATH) || {}
      yaml_content['defaults'] || {}
    else
      {}
    end
  end
  
  def self.camelize(str)
    str.split(/[-_]/).map(&:capitalize).join
  end
  
  def self.deep_merge(target, source)
    target.merge(source) do |key, oldval, newval|
      if oldval.is_a?(Hash) && newval.is_a?(Hash)
        deep_merge(oldval, newval)
      else
        newval
      end
    end
  end
  
  def self.build_config(key, pipeline_info)
    nf_core_name = pipeline_info['nf_core_name'] || key
    
    # Map category if possible
    category = pipeline_info['category']
    if category && CATEGORY_MAPPING[category]
      category = CATEGORY_MAPPING[category]
    end
    
    # Determine input type: manual config > auto-detect > default
    input_type = pipeline_info['input_type']
    if input_type.nil? || input_type == 'auto'
      # Auto-detect from schema (lazy - will be done when needed)
      input_type = nil  # Will be detected later
    end
    
    # Build config with convention defaults
    config = DEFAULTS.merge({
      'name' => "NfCore#{camelize(nf_core_name)}",
      'nf_core_name' => nf_core_name,
      'description' => pipeline_info['description'] || "nf-core/#{nf_core_name} pipeline",
      'analysis_category' => pipeline_info['analysis_category'] || category || 'nf-core',
      'default_version' => pipeline_info['default_version'] || pipeline_info['latest_version'] || 'master',
      'required_columns' => pipeline_info['required_columns'] || ['Name'],
      'input_type' => input_type,
      'custom_params' => pipeline_info['custom_params'] || [],
      'samplesheet_mapping' => pipeline_info['samplesheet_mapping'] || {}
    })
    
    # Merge any additional params from pipeline_info
    if pipeline_info['params']
      config['params'] = (config['params'] || {}).merge(pipeline_info['params'])
    end
    
    config
  end
  
  def self.generate_app_files
    load_pipelines.each do |key, pipeline_info|
      config = build_config(key, pipeline_info)
      
      class_name = "#{config['name']}App"
      file_path = File.join(LIB_DIR, "#{class_name}.rb")
      
      content = generate_class_code(class_name, config)
      File.write(file_path, content)
      puts "Generated #{file_path}"
    end
  end
  
  @@registered_classes = []
  
  def self.register_dynamic_apps
    puts "NfCoreAppFactory: register_dynamic_apps called"
    $stdout.flush
    pipelines = load_pipelines
    puts "NfCoreAppFactory: Found #{pipelines.size} pipelines from API"
    $stdout.flush
    
    newly_registered = 0
    pipelines.each do |key, pipeline_info|
      config = build_config(key, pipeline_info)
      class_name = "#{config['name']}App"
      
      # Always track the class name for nf-core apps
      unless @@registered_classes.include?(class_name)
        @@registered_classes << class_name
      end
      
      # Skip creating if class is already defined (static .rb file or previous registration)
      if Object.const_defined?(class_name)
        next
      end
      
      begin
        # Define class dynamically
        klass = create_dynamic_class(config)
        
        # Register as global constant
        Object.const_set(class_name, klass)
        newly_registered += 1
      rescue => e
        warn "NfCoreAppFactory: Failed to register #{class_name}: #{e.message}"
        # Remove from registered list if failed
        @@registered_classes.delete(class_name)
      end
    end
    
    puts "NfCoreAppFactory: Total nf-core apps tracked: #{@@registered_classes.size}, newly created: #{newly_registered}"
    $stdout.flush
  end
  
  def self.registered_class_names
    @@registered_classes
  end
  
  def self.create_dynamic_class(config)
    # Capture config in local variable for closure
    app_config = config.dup
    # Use fixed path accessible from all nodes (instead of ezRun package)
    r_app_path = '/srv/GT/analysis/masaomi/2026/FGCZ/nf_core_sushi_app_20260109/EzAppNfCoreGeneric.R'
    
    Class.new(SushiFabric::SushiApp) do
      # Need to include GlobalVariables to use run_RApp
      # But we need to override run method to call parent's run, not GlobalVariables' run
      include GlobalVariables
      
      # Override run to call SushiFabric::SushiApp's run method directly
      # (skipping GlobalVariables' run which requires an argument)
      define_method(:run) do
        SushiFabric::SushiApp.instance_method(:run).bind(self).call
      end
      
      define_method(:initialize) do
        super()
        @name = app_config['name']
        @analysis_category = app_config['analysis_category']
        @description = app_config['description']
        @required_columns = app_config['required_columns']
        @required_params = app_config['required_params'] || []
        @params['process_mode'] = 'DATASET'
        @params['nfcorePipeline'] = app_config['nf_core_name']
        @params['pipelineVersion'] = app_config['default_version']
        
        # Store input configuration for R script
        @params['inputType'] = app_config['input_type'] || 'samplesheet'
        
        # Store samplesheet mapping as JSON string for R
        if app_config['samplesheet_mapping'] && !app_config['samplesheet_mapping'].empty?
          @params['samplesheetMapping'] = app_config['samplesheet_mapping'].to_json
        end
        
        # Default params
        (app_config['params'] || {}).each do |k, v|
          @params[k] = v
        end
        
        @modules = ["Dev/jdk", "Tools/Nextflow"]
      end
      
      define_method(:next_dataset) do
        result_dir = File.join(@result_dir, "#{@params['name']}_result")
        
        dataset = {
          'Name' => @params['name'],
          'Result [File]' => result_dir,
          'MultiQC [Link]' => File.join(result_dir, 'multiqc', 'multiqc_report.html')
        }
        
        if @dataset && @dataset.first
          # Exclude input file columns (Read1, Read2 with any suffix like [File])
          # These should not be copied to output as they are input files
          inherit_cols = @dataset.first.keys.reject do |col|
            col == 'Name' || col =~ /^Read[12]/ || col == 'Species'
          end
          inherit_cols.each do |col|
            dataset[col] = @dataset.first[col]
          end
        end
        
        dataset
      end
      
      define_method(:grandchild_datasets) do
        []
      end
      
      define_method(:commands) do
        cmd = run_RApp('EzAppNfCoreGeneric')
        # Insert source() after the if block closes, before param = list()
        
        # Add Apptainer cache settings
        cache_settings = <<~SHELL
          export NXF_SINGULARITY_CACHEDIR=/misc/fgcz01/nextflow_apptainer_cache/
          export SINGULARITY_CACHEDIR=/misc/fgcz01/nextflow_apptainer_cache/
        SHELL
        
        cmd = cache_settings + cmd
        
        cmd.sub!("}\nparam = list()", "}\nsource('#{r_app_path}')\nparam = list()")
        cmd
      end
    end
  end
  
  def self.generate_class_code(class_name, config)
    # Use fixed path accessible from all nodes (instead of ezRun package)
    r_app_path = '/srv/GT/analysis/masaomi/2026/FGCZ/nf_core_sushi_app_20260109/EzAppNfCoreGeneric.R'
    
    <<~RUBY
      #!/usr/bin/env ruby
      # encoding: utf-8
      
      require 'sushi_fabric'
      require_relative 'global_variables'
      include GlobalVariables
      
      class #{class_name} < SushiFabric::SushiApp
        def initialize
          super
          @name = '#{config['name']}'
          @analysis_category = '#{config['analysis_category']}'
          @description =<<-EOS
      #{config['description']}
      EOS
          @required_columns = #{config['required_columns'].inspect}
          @required_params = #{config['required_params'] ? config['required_params'].inspect : '[]'}
          @params['process_mode'] = 'DATASET'
          @params['nfcorePipeline'] = '#{config['nf_core_name']}'
          @params['pipelineVersion'] = '#{config['default_version']}'
          
          # Default params
          #{(config['params'] || {}).map { |k, v| "@params['#{k}'] = #{v.inspect}" }.join("\n    ")}
          
          @modules = ["Dev/jdk", "Tools/Nextflow"]
        end
        
        def next_dataset
          result_dir = File.join(@result_dir, "\#{@params['name']}_result")
          
          dataset = {
            'Name' => @params['name'],
            'Result [File]' => result_dir,
            'MultiQC [Link]' => File.join(result_dir, 'multiqc', 'multiqc_report.html')
          }
          
          if @dataset && @dataset.first
            # Exclude input file columns (Read1, Read2 with any suffix like [File])
            # These should not be copied to output as they are input files
            inherit_cols = @dataset.first.keys.reject do |col|
              col == 'Name' || col =~ /^Read[12]/ || col == 'Species'
            end
            inherit_cols.each do |col|
              dataset[col] = @dataset.first[col]
            end
          end
          
          dataset
        end
        
        def grandchild_datasets
          []
        end
        
        def commands
          cmd = run_RApp('EzAppNfCoreGeneric')
          
          # Add Apptainer cache settings
          cache_settings = <<~SHELL
            export NXF_SINGULARITY_CACHEDIR=/misc/fgcz01/nextflow_apptainer_cache/
            export SINGULARITY_CACHEDIR=/misc/fgcz01/nextflow_apptainer_cache/
          SHELL
          
          cmd = cache_settings + cmd
          
          # Insert source() after the if block closes, before param = list()
          cmd.sub!("}\\nparam = list()", "}\\nsource('#{r_app_path}')\\nparam = list()")
          cmd
        end
      end
      
      if __FILE__ == $0
        usecase = #{class_name}.new
        usecase.project = "p1001"
        usecase.user = "sushi_lover"
        
        # Test run
        # usecase.test_run
      end
    RUBY
  end
end
