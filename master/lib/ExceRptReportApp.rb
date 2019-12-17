#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class ExceRptReportApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'ExceRptReport'
    @params['process_mode'] = 'DATASET'
    @params['name']
    @analysis_category = 'Count'
    @description =<<-EOS
    EOS
    ## modules
    @modules = ["Dev/R"]
  end
  def preprocess
    @random_string = (1..12).map{[*('a'..'z')].sample}.join
  end  
  def next_dataset
    report_file = File.join(@result_dir, "#{@params['name']}_ExceRptReport")
    report_link = File.join(report_file, '00index.html')
    dataset = {
      'Name'=>@dataset['Name'],
      'excerpt [File]'=>File.join(@result_dir, "#{@dataset['Name']}"),
      'Species'=>@dataset['Species'],
      'refBuild'=>@params['refBuild'],
    }
   dataset
  end
  def commands
    run_RApp("EzAppExceRptReport")
  end
end

if __FILE__ == $0
  run EzAppExceRptReport
  
end