#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class NfCoreAtacSeqApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'NfCoreAtacSeq'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'ATAC'
    @description =<<-EOS
    A ATAC-seq processing pipeline from NF-Core. <br/>
   <a href='https://nf-co.re/atacseq'>NF-Core ATAC-seq</a>
EOS
    @required_columns = ['Name','Read1','Read2','Species']
    @required_params = ['refBuild','readLength', 'peakStyle']
    @params['cores'] = '8'
    @params['ram'] = '100'
    @params['scratch'] = '200'
    @params['paired'] = true
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['readLength'] = ['150', '100', '50']
    @params['peakStyle'] = ['broad', 'narrow']
    @params['varStabilizationMethod'] = ['vst', 'rlogTransf'] 
    @params['name'] = 'NfCoreAtacSeq'
    @params['cmdOptions'] = ""
    @params['mail'] = ""
    @modules = ["Dev/jdk"]
  end
 def set_default_parameters
    if @params['paired']
      @required_columns<<  'Read2'
    end
  end
  def next_dataset
     report_file = File.join(@result_dir, "#{@params['name']}_results")
     report_link = File.join(@result_dir, "#{@params['name']}_results", "multiqc", "#{@params['peakStyle']}_peak", "multiqc_report.html")
     {'Name'=>@params['name'],
     'Species'=>(dataset = @dataset.first and dataset['Species']),
     'refBuild'=>@params['refBuild'],
     'Result [File]'=>report_file,
     'Report [Link]'=>report_link  
    }
  end
  def commands
    run_RApp('EzAppNfCoreAtacSeq')
  end
end
