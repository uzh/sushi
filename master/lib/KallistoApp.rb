#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class KallistoApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Kallisto'
    @analysis_category = 'Count'
    @description =<<-EOS
    <a href="https://pachterlab.github.io/kallisto/about">kallisto</a> is a program for quantifying abundances of transcripts from RNA-Seq data. It is based on the novel idea of pseudoalignment for rapidly determining the compatibility of reads with targets, without the need for alignment.
EOS
    @required_columns = ['Name','Read1','Species']
    @required_params = ['refBuild', 'paired', 'strandMode']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '100'
    @params['refBuild'] = ref_selector
    @params['paired'] = false
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['refFeatureFile'] = 'genes.gtf'
    @params['bootstrap-samples'] = 10
    @params['bootstrap-samples', 'description'] = 'number of bootstrap samples'
    @params['seed'] = 42
    @params['seed', 'description'] = 'seed for the bootstrap sampling'
    @params['fragment-length'] = 0
    @params['fragment-length', 'description'] = 'estimated average fragment length (required for single-end reads but should be set to 0 for paired-end reads)'
    @params['sd'] = 0
    @params['sd', 'description'] = 'estimated fragment length standard deviation (required for single-end reads but should be set to 0 for paired-end reads)'
    @params['bias'] = true
    @params['bias', 'description'] = 'perform sequence based bias correction'
    @params['pseudobam'] = false
    @params['pseudobam', 'description'] = 'generate a bam file with pseudoalignments'
    @params['transcriptFasta'] = ''
    @params['transcriptFasta', 'description'] = 'give full path of transcript fasta file; in that case the build is ignored; if it comes from trinity assembly the gene-isoform associations will be extracted and used'
    @params['transcriptTypes'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA', 'long_noncoding', 'short_noncoding', 'pseudogene']
    @params['transcriptTypes', 'multi_selection'] = true
    @params['trimAdapter'] = true
    @params['trimLeft'] = 0
    @params['trimRight'] = 0
    @params['minTailQuality'] = 0
    @params['minAvgQuality'] = 20
    @params['mail'] = ""
    @modules = ["Tools/samtools", "Aligner/kallisto", "QC/Flexbar", "QC/Trimmomatic", "Dev/R", "Tools/sambamba"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    else
      @required_params << 'fragment-length' << 'sd'
    end
  end
  def set_default_parameters
    @params['paired'] = dataset_has_column?('Read2')
    if dataset_has_column?('strandMode')
      @params['strandMode'] = @dataset[0]['strandMode']
    end
  end
  def next_dataset
    ds = {
      'Name'=>@dataset['Name'],
      'Count [File]'=>File.join(@result_dir, "#{@dataset['Name']}.txt"),
      'bootstrappedCount [File]'=>File.join(@result_dir, "#{@dataset['Name']}.h5"),
      'runInfo [File]'=>File.join(@result_dir, "#{@dataset['Name']}.json")
    }
    if @params['pseudobam']
      ds.merge!(
       'BAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam"),
       'BAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam.bai")
      )
    end
    ds.merge(
      'Species'=>@dataset['Species'],
      'refBuild'=>@params['refBuild'],
      'featureLevel'=>'isoform',
      'refFeatureFile'=>@params['refFeatureFile'],
      'strandMode'=>@params['strandMode'],
      'paired'=>@params['paired'],
      'Read Count'=>@dataset['Read Count'],
      'transcriptTypes'=>@params['transcriptTypes']
    ).merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppKallisto")
  end
end

if __FILE__ == $0
  run KallistoApp

end
