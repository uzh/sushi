#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class CellRangerARCApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'CellRangerARCCount'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
This wrapper runs <a href='https://support.10xgenomics.com/single-cell-ARC/software/pipelines/latest/using/count',>cellranger ARC count</a> in Single-library analysis mode.
    EOS
    @required_columns = ['Name','RawDataDir','Species','RawATACDataDir']
    @required_params = ['name', 'refBuild']
    @params['cores'] = '8'
    @params['ram'] = '60'
    @params['scratch'] = '200'
    @params['name'] = 'CellRangerARCCount'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['featureLevel'] = 'gene'
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'specify the commandline options for CellRanger ARC; do not specify any option that is already covered by the dedicated input fields'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Dev/R", "Aligner/CellRangerARC"]
    @inherit_tags = ["B-Fabric"]
  end
  def set_default_parameters
  end
  def next_dataset
    report_dir = File.join(@result_dir,"#{@dataset['Name']}")
    dataset = {
        'Name'=>@dataset['Name'],
        'Species'=>@dataset['Species'],
        'refBuild'=>@params['refBuild'],
        'refFeatureFile'=>@params['refFeatureFile'],
        'featureLevel'=>@params['featureLevel'],
        'ResultDir [File]'=>report_dir,
        'Report [Link]'=>File.join(report_dir, 'web_summary.html'),
        'CountMatrix [Link]'=>File.join(report_dir, 'filtered_peak_bc_matrix')
    }.merge(extract_columns(@inherit_tags))
    dataset
  end
  def commands
    run_RApp("EzAppCellRangerARC")
  end
end

if __FILE__ == $0

end
