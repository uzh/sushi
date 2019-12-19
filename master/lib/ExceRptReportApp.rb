#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class ExceRptReportApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'ExceRptReport'
    @required_columns = ['Name','excerpt', 'Species', 'refBuild']
    @params['process_mode'] = 'DATASET'
    @params['name'] = "Excerpt_Report"
    @analysis_category = 'QC'
    @description =<<-EOS
    EOS
    
    ## modules
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
    @modules = ["Dev/R"]
  end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@params['name'],
     'Species'=>(dataset = @dataset.first and dataset['Species']),
     'refBuild'=>@params['refBuild'],
     'Static Report [Link]'=>report_link,
     'Report [File]'=>report_file,
    }
  end
  def commands
    run_RApp("EzAppExceRptReport")
  end
end

if __FILE__ == $0
  run EzAppExceRptReport
  
end