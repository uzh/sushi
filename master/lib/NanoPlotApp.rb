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
    @required_columns = ['Name','Read1']
    @required_params = ['cores']
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '30'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '100'
    @params['scratch', "context"] = "slurm"
    @params['cmdOptions'] = ""
    @params['cmdOptions', "context"] = "NanoPlot"
    @params['mail'] = ""
    @modules = ["Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'NanoPlot Result [File]'=>File.join(@result_dir, "#{@dataset['Name']}"),
     'Report [Link]'=>File.join(@result_dir, "#{@dataset['Name']}", "#{@dataset['Name']}.NanoPlot-report.html")
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppNanoPlot",conda_env: "gi_nanoplot1.43.0")
  end
end

if __FILE__ == $0
  
end
