#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class NanoPlotApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'NanoPlot'
    @params['process_mode'] = 'SAMPLE'
    @analysis_category = 'QC'
    @description =<<-EOS
   NanoPlot: Plotting tool for long read sequencing data and alignments<br/>

    <a href='https://github.com/wdecoster/NanoPlot'/>Github web-site</a> 
EOS
    @required_columns = ['Name','Read1', 'Species']
    @required_params = ['name']
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '100'
    @params['name'] = 'NanoPlot_Result'
    @params['cmdOptions'] = ""
    @params['mail'] = ""
    @modules = ["Dev/R"]
  end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])    
    report_link = File.join(@result_dir, @params['name'], "#{@dataset['Name']}.NanoPlot-report.html")
    {'Name'=>@dataset['Name'],
     'Species'=>(dataset = @dataset.first and dataset['Species']),
     'Report [File]'=>report_file,
     'Report [Link]'=>report_link
    }
  end
  def commands
    run_RApp("EzAppNanoPlot",conda_env: "NanoPlot")
  end
end

if __FILE__ == $0
  
end
