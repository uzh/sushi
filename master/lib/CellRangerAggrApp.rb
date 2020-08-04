#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class CellRangerAggrApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'CellRangerAggr'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
This wrapper runs <a href='https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count',>cellranger aggr</a> in multi-library analysis mode.
EOS
    @required_columns = ['Name','CountMatrix']
    @required_params = ['name']
    @params['cores'] = '1'
    @params['ram'] = '4'
    @params['scratch'] = '100'
    @params['name'] = 'CellRangerAggr_Result'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['featureLevel'] = 'gene'
    @params['transcriptTypes'] = ''
    @params['transcriptTypes', 'multi_selection'] = true
    @params['transcriptTypes', 'selected'] = 0
    @params['normalize'] = ['mapped', 'none']
    @params['mail'] = ""
    @modules = ["Dev/R", "Aligner/CellRanger"]
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('refFeatureFile')
      @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    end
    @params['transcriptTypes'] = @dataset[0]['transcriptTypes'].to_s.split(',')
  end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    dataset = {
        'Name'=>@params['name'],
        'Condition'=>@dataset[0]['Condition'],
        'Species'=>@dataset[0]['Species'],
        'refBuild'=>@params['refBuild'],
        'refFeatureFile'=>@params['refFeatureFile'],
        'featureLevel'=>@params['featureLevel'],
        'transcriptTypes'=>@params['transcriptTypes'],
        'ResultDir [File]'=>report_file,
        'Report [Link]'=>File.join(report_file, 'web_summary.html'),
        'CountMatrix [Link]'=>File.join(report_file, 'filtered_feature_bc_matrix')
      }
  end
  def commands
    run_RApp("EzAppCellRangerAggr")
  end
end

