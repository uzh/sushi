#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class CellRangerMultiApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'CellRangerMulti'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
This wrapper runs <a href='https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count',>cellranger count</a> in Single-library analysis mode.
    EOS
    @required_columns = ['Name','RawDataDir','Species']
    @required_params = ['name', 'refBuild']
    @params['cores'] = ['8', '12', '16']
    @params['ram'] = ['60', '40', '80']
    @params['scratch'] = '200'
    @params['name'] = 'CellRangerMulti'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['featureLevel'] = 'gene'
    @params['TenXLibrary'] = ['GEX', 'VDJ', 'FeatureBarcoding', 'Multiplexing']
    @params['TenXLibrary', 'description'] = 'Which 10X libraries? Note: Not all library types can be processed simultaneously. For GEX or GEX+FeatureBarcoding, the CellRangerCount app is recommended'
    @params['TenXLibrary', 'multi_selection'] = true
    @params['TenXLibrary', 'selected'] = ['GEX'] 
    @params['FeatureBarcodeFile'] = ''
    @params['FeatureBarcodeFile', 'file_upload'] = true
    @params['FeatureBarcodeFile', 'description'] = '(e.g. for CITEseq)'
    @params['MultiplexBarcodeFile'] = {'select'=>''}
    @params['MultiplexBarcodeFile'] =  {'select'=>''}
    Dir["/srv/GT/databases/10x/CMO_files/*"].sort.select{|design| File.file?(design)}.each do |dir|
      @params['MultiplexBarcodeFile'][File.basename(dir)] = File.basename(dir)
    end
    @params['MultiplexBarcodeFile', 'description'] = 'Used when CellPlex libraries. New files needs to be installed under /srv/GT/databases/10x/CMO_files'
    @params['SampleMultiplexBarcodeFile'] = ''
    @params['SampleMultiplexBarcodeFile', 'file_upload'] = true
    @params['SampleMultiplexBarcodeFile', 'description'] = 'For assigning multiplexing barcode IDs to samples'
    @params['includeIntrons'] = true
    @params['includeIntrons', 'description'] = 'set to false to reproduce the default behavior in cell ranger v6 and earlier'
    @params['expectedCells'] = ''
    @params['expectedCells', 'description'] = 'Expected number of recovered cells. Leave this free to let cellranger estimate the cell number.'
    @params['transcriptTypes'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA', 'long_noncoding', 'short_noncoding', 'pseudogene']
    @params['transcriptTypes', 'multi_selection'] = true
    @params['transcriptTypes', 'selected'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA']
    @params['controlSeqs'] = ''
    @params['controlSeqs', 'description'] = 'The extra control sequences (such as spikein sequences) available in https://fgcz-gstore.uzh.ch/reference/controlSeqs.fa'
    @params['bamStats'] = true
    @params['bamStats', 'description'] = 'Compute stats per cell from the bam file?'
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'specify the commandline options for CellRanger (e.g. --include-introns for single nuclei data); do not specify any option that is already covered by the dedicated input fields'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Tools/seqtk", "Dev/R/4.3.0", "Dev/Python/3.8.3", "Tools/samtools"]
    @params['CellRangerVersion'] = ["Aligner/CellRanger/7.1.0"]
    @inherit_tags = ["Factor", "B-Fabric"]
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
      'transcriptTypes'=>@params['transcriptTypes'],
      'ResultDir [File]'=>report_dir,
      'Report [Link]'=>File.join(report_dir, 'web_summary.html'),
      'CountMatrix [Link]'=>File.join(report_dir, 'filtered_feature_bc_matrix'),
      'Read Count'=>@dataset['Read Count']
    }.merge(extract_columns(@inherit_tags))
    dataset
  end
  def commands
    command = "module load  #{@params["CellRangerVersion"]}\n"
    command << run_RApp("EzAppCellRangerMulti")
  end
end

if __FILE__ == $0

end
