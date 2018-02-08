#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class AtaqvApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'Ataqv'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'ATAC'
    @description =<<-EOS
    A toolkit for QC and visualization of ATAC-seq results. <br/>

    <a href='https://github.com/ParkerLab/ataqv'/>Ataqv web-site</a>
EOS
    @required_columns = ['Name','BAM', 'Species', 'refBuild']
    @required_params = ['name']
    @params['cores'] = '8'
    @params['ram'] = '16'
    @params['scratch'] = '50'
    @params['refBuild'] = ref_selector
    @params['paired'] = true
    @params['name'] = 'Ataqv'
    @params['mail'] = ""
    @modules = ["Dev/R", "Tools/ataqv"]
  end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, "index.html")
    {'Name'=>@params['name'],
     'Report [File]'=>report_file,
     'Html [Link]'=>report_link,
    }
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
  end
  def commands
    run_RApp("EzAppAtaqv")
  end
end

if __FILE__ == $0
  
end
