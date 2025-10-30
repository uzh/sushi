#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class EzPyzAnndataReportApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'EzPyzAnndataReportApp'
    @params['process_mode'] = 'SAMPLE'
    @analysis_category = 'Spatial'
    @description =<<-EOS
    A Anndata report app for Spatial Trancriptomics data.<br/>
    This app is running fully on Python, therefore, check the implementation in ezPyz: https://github.com/fgcz/EzPyzApps . There you can find the required information and the implementation of the app.<br/>
    The compiled code lives in the conda environment `tmp_anndata_report`<br/>
    EOS
    @required_columns = ['Name','Anndata']
    @required_params = ['name']
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '30'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '100'
    @params['scratch', "context"] = "slurm"
    @params['name'] = 'AnndataReportApp'
    @params['mail'] = ""
    @modules = []
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def next_dataset    
    dir_name = "#{@params['name']}_#{@dataset['Name']}"
    report_dir = File.join(@result_dir, dir_name)
    {'Name'=> dir_name,
    'Report Folder [File]'=>report_dir,
    'Report [Link]'=>File.join(report_dir, 'report/AnnData_quality_report.html')
    }
  end
  def commands
    run_PyApp("AnndataReport",conda_env: 'tmp_anndata_report')
  end
end
