#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class NfCoreRnaseqApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'NfCoreRnaseq'
    @analysis_category = 'Transcriptomics'
    @description =<<-EOS
RNA sequencing analysis pipeline using STAR, RSEM, HISAT2 or Salmon with gene/isoform counts and extensive quality control.
<a href='https://api.github.com/repos/nf-core/rnaseq'>nf-core/rnaseq</a>

EOS
    @required_columns = ["Name", "Read1"]
    @required_params = []
    @params['process_mode'] = 'DATASET'
    @params['nfcorePipeline'] = 'rnaseq'
    @params['pipelineVersion'] = '3.22.2'
    
    # Default params
    @params['cores'] = 8
    @params['ram'] = 30
    @params['scratch'] = 100
    
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
    
    # Add Apptainer cache settings
    cache_settings = <<~SHELL
      export NXF_SINGULARITY_CACHEDIR=/misc/fgcz01/nextflow_apptainer_cache/
      export SINGULARITY_CACHEDIR=/misc/fgcz01/nextflow_apptainer_cache/
    SHELL
    
    cmd = cache_settings + cmd
    
    # Insert source() after the if block closes, before param = list()
    r_app_path = '/srv/GT/analysis/masaomi/2026/FGCZ/nf_core_sushi_app_20260109/EzAppNfCoreGeneric.R'
    cmd.sub!("}\nparam = list()", "}\nsource('#{r_app_path}')\nparam = list()")
    cmd
  end
end

if __FILE__ == $0
  usecase = NfCoreRnaseqApp.new
  usecase.project = "p1001"
  usecase.user = "sushi_lover"
  
  # Test run
  # usecase.test_run
end
