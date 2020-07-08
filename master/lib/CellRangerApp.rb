#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class CellRangerApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'CellRangerCount'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
This wrapper runs <a href='https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count',>cellranger count</a> in Single-library analysis mode.
    EOS
    @required_columns = ['Name','RawDataDir','Species']
    @required_params = ['name', 'refBuild']
    @params['cores'] = '8'
    @params['ram'] = '62'
    @params['scratch'] = '200'
    @params['name'] = 'CellRangerCount'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['featureLevel'] = 'gene'
    @params['TenXLibrary'] = ['GEX', 'VDJ']
    @params['TenXLibrary', 'description'] = 'Which 10X library? GEX or VDJ.'
    @params['scMode'] = ['SC', 'SN']
    @params['scMode', 'description'] = 'Single-cell or single-nuclei?'
    @params['chemistry'] = ['auto', 'threeprime', 'fiveprime', 'SC3Pv1', 'SC3Pv2', 'SC3Pv3', 'SC5P-PE', 'SC5P-R2']
    @params['chemistry', 'description'] = 'Assay configuration. NOTE: by default the assay configuration is detected automatically, which is the recommended mode.'
    @params['transcriptTypes'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA', 'long_noncoding', 'short_noncoding', 'pseudogene']
    @params['transcriptTypes', 'multi_selection'] = true
    @params['transcriptTypes', 'selected'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA']
    @params['controlSeqs'] = ''
    @params['controlSeqs', 'description'] = 'The extra control sequences (such as spikein sequences) available in https://fgcz-gstore.uzh.ch/reference/controlSeqs.fa'
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'specify the commandline options for CellRanger; do not specify any option that is already covered by the dedicated input fields'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Dev/R", "Aligner/CellRanger"]
    @inherit_tags = ["B-Fabric"]
  end
  def set_default_parameters
  end
  def next_dataset
    report_dir = File.join(@result_dir,"#{@dataset['Name']}")
    if @params['TenXLibrary'] == "VDJ"
      dataset = {
        'Name'=>@dataset['Name'],
        'Species'=>@dataset['Species'],
        'refBuild'=>@params['refBuild'],
        'refFeatureFile'=>@params['refFeatureFile'],
        'featureLevel'=>@params['featureLevel'],
        'ResultDir [File]'=>report_dir,
        'Report [Link]'=>File.join(report_dir, 'web_summary.html')
      }.merge(extract_columns(@inherit_tags))
    else
      dataset = {
        'Name'=>@dataset['Name'],
        'Condition'=>@dataset['Condition'],
        'Species'=>@dataset['Species'],
        'refBuild'=>@params['refBuild'],
        'refFeatureFile'=>@params['refFeatureFile'],
        'featureLevel'=>@params['featureLevel'],
        'transcriptTypes'=>@params['transcriptTypes'],
        'ResultDir [File]'=>report_dir,
        'Report [Link]'=>File.join(report_dir, 'web_summary.html'),
        'CountMatrix [Link]'=>File.join(report_dir, 'filtered_feature_bc_matrix')
      }.merge(extract_columns(@inherit_tags))
    end
    dataset
  end
  def commands
    run_RApp("EzAppCellRanger")
  end
end

if __FILE__ == $0

end
