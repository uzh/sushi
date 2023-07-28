#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class BDRhapsodyApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'BDRhapsodySequenceAnalysis'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
      This wrapper runs a <a href='https://scomix.bd.com/hc/en-us/categories/360000838932-Resource-Library',>CWL workflow</a> for the analysis of BD Single-Cell Multiomics.
    EOS
    @required_columns = ['Name', 'Read1', 'Read2', 'Species']
    @required_params = ['name', 'refBuild']
    @params['cores'] = ['8', '12', '16']
    @params['ram'] = ['60', '40', '80']
    @params['scratch'] = '200'
    @params['name'] = 'BDRhapsodySequenceAnalysis'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['featureLevel'] = 'gene'
    @params['abSeqReference'] = ''
    @params['abSeqReference', 'file_upload'] = true
    @params['abSeqReference', 'description'] = 'FASTA AbSeq reference file generated from https://abseq-ref-gen.genomics.bd.com/'
    @params['enableRefinedPutativeCellCalling'] = false
    @params['enableRefinedPutativeCellCalling', 'description'] = 'Enable use of refined putative cell calling algorithm for cell calling if set to True. By default, putative cells are determined using only the basic algorithm (minimum second derivative along the cumulative reads curve).'
    @params['exactCellCount'] = ''
    @params['exactCellCount', 'description'] = 'Set a specific number (>=1) of cells as putative, based on those with the highest error-corrected read count.'
    @params['expectedCellCount'] = ''
    @params['expectedCellCount', 'description'] = 'Guide the basic putative cell calling algorithm by providing an estimate of the number of cells expected. Usually this can be the number of cells loaded into the Rhapsody cartridge.'
    @params['generateBamOutput'] = false
    @params['excludeIntronicReads'] = false
    @params['putativeCellCalling'] = ['mRNA', 'AbSeq']
    @params['putativeCellCalling', 'description'] = 'Specify the data to be used for putative cell calling'
    @params['targetedReference'] = ''
    @params['targetedReference', 'file_upload'] = true
    @params['targetedReference', 'description'] = 'This is an mRNA reference file in fasta format. This is a pre-designed, supplemental, or custom panel.'
    @params['transcriptTypes'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA', 'long_noncoding', 'short_noncoding', 'pseudogene']
    @params['transcriptTypes', 'multi_selection'] = true
    @params['transcriptTypes', 'selected'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA']
    @params['sampleTagsVersion'] = ['', 'Single-Cell Multiplex Kit - Human', 'Single-Cell Multiplex Kit - Mouse', 'Single-Cell Multiplex Kit - Flex']
    @params['tagNames'] = ''
    @params['tagNames', 'description'] = "For a multiplexed samples run only. Use the format '[Sample Tag number]-[sample name]'. Comma-separation for multiple samples. No spaces or forward slashes allowed. Example: 3-Foo,4-Bar"
    @params['vdjSpeciesVersion'] = ['', 'Human VDJ - BCR and TCR', 'Human VDJ - BCR only', 'Human VDJ - TCR only', 'Mouse VDJ - BCR and TCR', 'Mouse VDJ - BCR only', 'Mouse VDJ - TCR only']
    @params['controlSeqs'] = ''
    @params['controlSeqs', 'description'] = 'The extra control sequences (such as spikein sequences) available in https://fgcz-gstore.uzh.ch/reference/controlSeqs.fa'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Dev/R/4.3.0"]
    @inherit_tags = ["Order Id", "Factor", "B-Fabric"]
  end
  def set_default_parameters
  end
  def next_dataset
    report_dir = File.join(@result_dir, "#{@dataset['Name']}")
    dataset = {
      'Name'=>@dataset['Name'],
      'Species'=>@dataset['Species'],
      'refBuild'=>@params['refBuild'],
      'refFeatureFile'=>@params['refFeatureFile'],
      'featureLevel'=>@params['featureLevel'],
      'transcriptTypes'=>@params['transcriptTypes'],
      'ResultDir [File]'=>report_dir,
      'Report [Link]'=>File.join(report_dir, "#{@dataset['Name']}_Pipeline_Report.html"),
      'CountMatrix [Link]'=>File.join(report_dir, "#{@dataset['Name']}_DBEC_MolsPerCell_MEX"),
      'Read Count'=>@dataset['Read Count']
    }.merge(extract_columns(@inherit_tags))
    dataset
  end
  def commands
    run_RApp("EzAppBDRhapsody", conda_env: 'seven-bridges')
  end
end

if __FILE__ == $0

end
