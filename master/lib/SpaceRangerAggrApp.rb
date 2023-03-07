#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class SpaceRangerAggrApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'SpaceRangerAggr'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'SpatialTrx'
    @description =<<-EOS
This wrapper runs <a href='https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/using/aggregate',>spaceranger aggr</a> in multi-library analysis mode.
EOS
    @required_columns = ['Name','CountMatrix']
    @required_params = ['name']
    @params['cores'] = '1'
    @params['ram'] = ['7', '15']
    @params['scratch'] = ['100', '150']
    @params['name'] = 'SpaceRangerAggr_Result'
    @params['normalize'] = ['mapped', 'none']
    @params['mail'] = ""
    @modules = ["Dev/R", "Aligner/SpaceRanger"]
  end
  def set_default_parameters
  end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    dataset = {
        'Name'=>File.join(@dataset[0]['Name'],@dataset[1]['Name']),
        'Species'=>@dataset[0]['Species'],
        'refBuild'=>@dataset[0]['refBuild'],
        'refFeatureFile'=>@dataset[0]['refFeatureFile'],
        'featureLevel'=>@dataset[0]['featureLevel'],
        'transcriptTypes'=>@dataset[0]['transcriptTypes'],
        'ResultDir [File]'=>report_file,
        'Report [Link]'=>File.join(report_file, 'web_summary.html'),
        'CountMatrix [Link]'=>File.join(report_file, 'count', 'filtered_feature_bc_matrix')
      }
  end
  def commands
    run_RApp("EzAppSpaceRangerAggr")
  end
end

