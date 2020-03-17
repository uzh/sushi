#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class CellRangerATACAggrApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'CellRangerATACAggr'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
This wrapper runs <a href='https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/using/count',>cellranger atac aggr</a> in multi-library analysis mode.
EOS
    @required_columns = ['Name','CountMatrix']
    @required_params = ['name']
    @params['cores'] = '1'
    @params['ram'] = '4'
    @params['scratch'] = '100'
    @params['name'] = 'CellRangerATACAggr_Result'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['featureLevel'] = 'gene'
    @params['normalize'] = ['depth', 'signal', 'none']
    @params['mail'] = ""
    @modules = ["Dev/R"]
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('refFeatureFile')
      @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    end
  end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    {'Name'=>@params['name'],
     'refBuild'=>@params['refBuild'],
     'refFeatureFile'=>@params['refFeatureFile'],
     'featureLevel'=>@params['featureLevel'],
     'ResultDir [File]'=>report_file,
     'Report [Link]'=>File.join(report_file, 'web_summary.html')
    }
  end
  def commands
    run_RApp("EzAppCellRangerATACAggr")
  end
end
