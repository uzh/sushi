#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class VPipeApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'VPipe'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Variants'
    @description =<<-EOS
    V-pipe: A bioinformatics pipeline for viral sequencing data <br/>
    This app supports currently only SARS-CoV2 data <br/>

    <a href='https://github.com/cbg-ethz/V-pipe'/>Github web-site</a>
    
EOS
    @required_columns = ['Name','Read1','Read2','Species']
    @required_params = ['name', 'paired']
    @params['cores'] = '24'
    @params['ram'] = '50'
    @params['scratch'] = '100'
    @params['paired'] = true
    @params['name'] = 'VPIPE_Result'
    @params['readLength'] = '150'
    @params['samplePrefix'] = '210531_SARS_BAG_order25068_'
    @params['samplePrefix', 'description'] = 'prefix to remove from sample name'
    @params['cmdOptions'] = ""
    @params['mail'] = ""
    @modules = ["Dev/R", "Tools/samtools", "Tools/BEDTools"]
  end
 def set_default_parameters
    @params['paired'] = dataset_has_column?('Read2')
  end
  def preprocess
    if @params['paired']
      @required_columns<<  'Read2'
    end
  end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])    
    {'Name'=>@params['name'],
     'Report [File]'=>report_file,
    }
  end
  def commands
    run_RApp("EzAppVPipe",conda_env: "V-pipe")
  end
end

if __FILE__ == $0
  
end
