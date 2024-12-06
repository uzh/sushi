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
    @params['transcriptTypes', 'selected'] = ['protein_coding']
    @params['transcriptTypes', 'multi_selection'] = true

    # trimming options
    # general
    @params['trimAdapter', 'hr'] = true
    @params['trimAdapter'] = true
    # Fastp
    ## trimming
    @params['trim_front1'] = '0'
    @params['trim_front1','description'] = 'trimming how many bases in front for read1 (and read2), default is 0.'
    @params['trim_tail1'] = '0'
    @params['trim_tail1','description'] = 'trimming how many bases in tail for read1 (and read2), default is 0.'
    @params['cut_front'] = false
    @params['cut_front','description'] = 'move a sliding window from front (5p) to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.'
    @params['cut_front_window_size'] = '4'
    @params['cut_front_mean_quality'] = '20'
    @params['cut_tail'] = false
    @params['cut_tail','description'] = 'move a sliding window from tail (3p) to front, drop the bases in the window if its mean quality < threshold, stop otherwise.'
    @params['cut_tail_window_size'] = '4'
    @params['cut_tail_mean_quality'] = '20'
    @params['cut_right'] = false
    @params['cut_right','description'] = 'move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, and then stop.'
    @params['cut_right_window_size'] = '4'
    @params['cut_right_mean_quality'] = '20'    
    @params['average_qual'] = '0'
    @params['average_qual','description'] = 'if one read average quality score <avg_qual>, then this read/pair is discarded. Default 0 means no requirement'
    @params['max_len1'] = '0'
    @params['max_len1','description'] = 'if read1 is longer than max_len1, then trim read1 at its tail to make it as long as max_len1. Default 0 means no limitation. If two reads are present, the same will apply to read2.'
    @params['max_len2'] = '0'
    @params['max_len2','description'] = 'if read1 is longer than max_len2, then trim read2 at its tail to make it as long as max_len1. Default 0 means no limitation.'
    @params['poly_x_min_len'] = '10'
    @params['poly_x_min_len','description'] = 'the minimum length to detect polyX in the read tail. 10 by default.'
    @params['length_required'] = '18'
    @params['length_required','description'] = 'reads shorter than length_required will be discarded.'
    @params['cmdOptionsFastp'] = ''
    ## additional commands
    @params['markDuplicates'] = false
    @params['markDuplicates', 'description'] = 'should duplicates be marked with picard. It is recommended for ChIP-seq and ATAC-seq data.'
      
    @params['mail'] = ""
    @modules = ["Tools/samtools", "Aligner/kallisto/0.46.1", "QC/fastp", "Dev/R"]
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
      'runInfo [File]'=>File.join(@result_dir, "#{@dataset['Name']}.json"),
      'PreprocessingLog [File]'=>File.join(@result_dir, "#{@dataset['Name']}_preprocessing.log"),
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
