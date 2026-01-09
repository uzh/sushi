require 'yaml'
require 'fileutils'
require_relative 'nf_core_info_fetcher'

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
  
  def self.generate_class_code(class_name, config)
    # R script path logic to be embedded in the generated file
    # We want the path to be resolved at runtime relative to the App file
    r_app_path_code = "'\#{File.expand_path('../R/EzAppNfCoreGeneric.R', __FILE__)}'"
    
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
            inherit_cols = @dataset.first.keys - ['Name', 'Read1', 'Read2', 'Species']
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
          cmd.sub!("}\\nparam = list()", "}\\nsource(#{r_app_path_code})\\nparam = list()")
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
