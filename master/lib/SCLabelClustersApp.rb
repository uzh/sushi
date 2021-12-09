#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

require 'csv'
class SCLabelClustersApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'SCLabelClusters'
    @description =<<-EOS
Updated clustering results from SCOneSample application: new annotations and eventual merging included.
EOS
    @analysis_category = 'SingleCell'
    @params['cores'] = '4'
    @params['ram'] = '30'
    @params['scratch'] = '50'

    @params['clusterInfo'] = ''
    @params['clusterInfo', 'file_upload'] = true
    @required_columns = ['Name', 'SC Cluster Report', 'SC H5']
    #@required_columns = ['Name']
    @required_params = ['clusterInfo']
    #@required_params = []
    #@inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
    @inherit_tags = []
  end
  def next_dataset
    report_file_CellType = File.join(@result_dir, 'SC_CellType_Report')
    report_file_H5= File.join(@result_dir, 'SC_H5_Report')
    {
      'Name'=>@dataset['Name'],
      'SC CellType Report [File]'=>report_file_CellType,
      'SC H5 [File]'=>report_file_H5
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppSCLabelClusters")
  end
end

