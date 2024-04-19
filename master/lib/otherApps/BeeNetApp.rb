#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class BeeNetApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'BeeNet'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
      This wrapper runs <a href='https://honeycombbio.zendesk.com/hc/en-us/sections/4408705359643-BeeNet',>BeeNet</a>, a custom software designed to process data from paired-end Illumina® sequencing of single-cell RNA-seq libraries produced by the HIVE scRNAseq Processing Kit.
    EOS
    @required_columns = ['Name', 'Read1', 'Read2', 'Species']
    @required_params = ['name', 'refBuild', 'numBarcodes']
    @params['cores'] = ['8', '12', '16']
    @params['ram'] = ['60', '40', '80']
    @params['scratch'] = '200'
    @params['name'] = 'BeeNet'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['featureLevel'] = 'gene'
    @params['transcriptTypes'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA', 'long_noncoding', 'short_noncoding', 'pseudogene']
    @params['transcriptTypes', 'multi_selection'] = true
    @params['transcriptTypes', 'selected'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA']
    @params['numBarcodes'] = '5000'
    @params['numBarcodes', 'description'] = 'This is the expected number of barcodes in a given sample. Count matrix generation uses this, after sorting by count of reads mapped to genes, as the number of cell barcodes in the sample. Must be an integer. We recommend this number be ~40% of the starting input cell number.'
    @params['controlSeqs'] = ''
    @params['controlSeqs', 'description'] = 'The extra control sequences (such as spikein sequences) available in https://fgcz-gstore.uzh.ch/reference/controlSeqs.fa'
    @params['specialOptions'] = ''
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'specify additional command-line options for beenet analyze (e.g. --qual-bc for changing the maximum allowed mismatches in the barcode); do not specify any option that is already covered by the dedicated input fields'
    @params['BeeNetVersion'] = ['1.1.2']
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
      'CountMatrix [Link]'=>File.join(report_dir, "TCM_counts_MEX"),
      'Read Count'=>@dataset['Read Count']
    }.merge(extract_columns(@inherit_tags))
    dataset
  end
  def commands
    command = "module load Aligner/BeeNet/#{@params["BeeNetVersion"]}\n"
    command << run_RApp("EzAppBeeNet")
  end
end

if __FILE__ == $0

end
