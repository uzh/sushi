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
    @required_columns = ['Name', 'Species', 'refBuild', 'SC Cluster Report', 'SC H5']
    @required_params = ['clusterInfo']
    #@inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
    @inherit_tags = []
  end
  def next_dataset
    report_file = File.join(@result_dir, "#{@dataset['Name']}_Celltype_Report")
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@dataset['Name'],
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'refFeatureFile'=>@params['refFeatureFile'],
     'Static Report [Link]'=>report_link,
     'SC Celltype Report [File]'=>report_file,
     'SC H5 [File]'=>File.join(report_file, "sce_h5"),
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppSCLabelClusters")
  end
end

