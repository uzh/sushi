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
    A Anndata report app for Spatial Trancriptomics data<br/>
    EOS
    @required_columns = ['Name','Anndata']
    @required_params = ['name']
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '100'
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
