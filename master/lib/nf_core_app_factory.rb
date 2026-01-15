require 'yaml'
require 'fileutils'
require_relative 'nf_core_info_fetcher'
require_relative 'global_variables'

module NfCoreAppFactory
  PIPELINE_CONFIG_PATH = File.expand_path('../../config/nf_core_pipelines.yml', __FILE__)
  LIB_DIR = File.expand_path('..', __FILE__)
  
  DEFAULTS = {
    'analysis_category' => 'Misc',
    'required_columns' => ['Name', 'Read1'],
    'params' => { 'cores' => 8, 'ram' => 30, 'scratch' => 100 }
  }
  
  CATEGORY_MAPPING = {
    'rna-seq' => 'Transcriptomics',
    'dna-seq' => 'Genomics',
    'chip-seq' => 'Epigenetics',
    'metagenomics' => 'Metagenomics',
    'quality-control' => 'QC'
  }
  
  def self.load_pipelines
    return {} unless File.exist?(PIPELINE_CONFIG_PATH)
    YAML.load_file(PIPELINE_CONFIG_PATH)['pipelines'] || {}
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
  
  def self.build_config(key, user_config)
    nf_core_name = user_config['nf_core_name'] || key
    
    # Fetch missing info from nf-core APIs
    fetched = NfCoreInfoFetcher.fetch_all(nf_core_name)
    
    # Map category if possible
    category = fetched[:category]
    if category && CATEGORY_MAPPING[category]
      category = CATEGORY_MAPPING[category]
    end
    
    # Build config with convention defaults
    config = DEFAULTS.merge({
      'name' => "NfCore#{camelize(nf_core_name)}",
      'nf_core_name' => nf_core_name,
      'description' => fetched[:description],
      'analysis_category' => category || 'Misc',
      'default_version' => fetched[:latest_version],
      'required_columns' => fetched[:required_columns]
    })
    
    # User config overrides everything
    deep_merge(config, user_config)
  end
  
  def self.generate_app_files
    load_pipelines.each do |key, user_config|
      config = build_config(key, user_config)
      
      class_name = "#{config['name']}App"
      file_path = File.join(LIB_DIR, "#{class_name}.rb")
      
      content = generate_class_code(class_name, config)
      File.write(file_path, content)
      puts "Generated #{file_path}"
    end
  end
  
  @@registered_classes = []
  
  def self.register_dynamic_apps
    load_pipelines.each do |key, user_config|
      config = build_config(key, user_config)
      class_name = "#{config['name']}App"
      
      # Skip if class is already defined
      next if Object.const_defined?(class_name)
      
      # Define class dynamically
      klass = create_dynamic_class(config)
      
      # Register as global constant
      Object.const_set(class_name, klass)
      @@registered_classes << class_name
    end
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
