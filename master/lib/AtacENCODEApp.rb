#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class AtacENCODEApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'AtacENCODE'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'ATAC'
    @description =<<-EOS
    A ATAC-seq and DNase-seq processing pipeline from ENCODE <br/>
    <a href='https://github.com/kundajelab/atac_dnase_pipelines'/>Github web-site</a>
    When the job fails, please go to Jobs to resubmit the job. 
    It will make sure the jobs reruns from the failed point.
EOS
    @required_columns = ['Name','Read1','Read2']
    @required_params = ['name', 'paired','Species']
    @params['cores'] = '8'
    @params['ram'] = '16'
    @params['scratch'] = '100'
    @params['paired'] = true
    @params['name'] = 'AtacENCODE'
    @params['cmdOptions'] = ""
    @params['mail'] = ""
    @modules = []
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
    report_link = File.join(report_file, "out", "#{@params['name']}_report.html")
    {'Name'=>@params['name'],
     'Report [File]'=>report_file,
     'Html [Link]'=>report_link,
    }
  end
  def commands
    run_RApp("EzAppAtacENCODE")
  end
end

if __FILE__ == $0
  
end
