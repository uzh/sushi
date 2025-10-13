#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class XeniumQCApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'XeniumQC'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Spatial'
    @description =<<-EOS
Multi-Sample Quality Control for Xenium Spatial Transcriptomics<br/>
Aggregates and visualizes QC metrics from multiple Xenium runs/regions.<br/>
Each region is treated as a separate sample in the report.
    EOS
    @required_columns = ['Name','XeniumPath']
    @required_params = ['name']
    @params['cores'] = '1'
    @params['ram'] = '30'
    @params['scratch'] = '100'
    @params['name'] = 'XeniumQC'
    @params['mail'] = ""
    @modules = ["Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def next_dataset
    report_dir = File.join(@result_dir, @params['name'])
    {'Name'=>@params['name'],
     'ReportData [File]'=>report_dir,
     'Report [Link]'=>File.join(report_dir, '00index.html'),
     'Species'=>(dataset = @dataset.first and dataset['Species'])
    }
  end
  def commands
    run_RApp("EzAppXeniumQC")
  end
end
