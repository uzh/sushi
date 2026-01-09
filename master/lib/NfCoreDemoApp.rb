#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class NfCoreDemoApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'NfCoreDemo'
    @analysis_category = 'QC'
    @description =<<-EOS
Demo pipeline from nf-core for testing.
<a href='https://nf-co.re/demo'>nf-core/demo</a>

EOS
    @required_columns = ["Name", "Read1"]
    @required_params = []
    @params['process_mode'] = 'DATASET'
    @params['nfcorePipeline'] = 'demo'
    @params['pipelineVersion'] = '1.0.2'
    
    # Default params
    @params['cores'] = 4
    @params['ram'] = 16
    @params['scratch'] = 50
    
    @modules = ["Dev/jdk", "Tools/Nextflow"]
  end
  
  def next_dataset
    result_dir = File.join(@result_dir, "#{@params['name']}_result")
    
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
    cmd.sub!("}\nparam = list()", "}\nsource('#{File.expand_path('../R/EzAppNfCoreGeneric.R', __FILE__)}')\nparam = list()")
    cmd
  end
end

if __FILE__ == $0
  usecase = NfCoreDemoApp.new
  usecase.project = "p1001"
  usecase.user = "sushi_lover"
  
  # Test run
  # usecase.test_run
end
