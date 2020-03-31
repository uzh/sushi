#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class SalmonApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Salmon'
    @analysis_category = 'Count'
    @description =<<-EOS
	 <a href="https://github.com/COMBINE-lab/salmon">Salmon</a> is a wicked-fast program to produce a highly-accurate, transcript-level quantification estimates from RNA-seq data.
EOS
    @required_columns = ['Name','Read1','Species']
    @required_params = ['paired', 'strandMode']
    @params['paired'] = false
    @params['strandMode'] = ['both', 'sense', 'antisense']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '100'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['fldMean'] = 0
    @params['fldMean', 'description'] = 'expected mean fragment length for single-end reads'
    @params['fldSD'] = 0
    @params['fldSD', 'description'] = 'expected standard deviation of the fragment length distribution for single-end reads'
    @params['useVBOpt'] = false
    @params['useVBOpt', 'description'] = 'use the variational Bayesian EM algorithm rather than the "standard" EM algorithm to optimize abundance estimates'
    @params['numBootstraps'] = 0
    @params['numBootstraps', 'description'] = 'number of bootstrapped samples; by default, Salmon does not use bootstrapping'
    @params['numGibbsSamples'] = 0
    @params['numGibbsSamples', 'description'] = 'an alternative to bootstrapping for estimating the variance in abundance estimates'
    @params['seqBias'] = false
    @params['seqBias', 'description'] = 'enable sequence bias correction'
    @params['gcBias'] = false
    @params['gcBias', 'description'] = 'enable GC bias correction'
    @params['posBias'] = false
    @params['posBias', 'description'] = 'enable position bias correction'
    @params['specialParams'] = ''
    @params['specialParams', 'description'] = 'additional command line parameters to pass to Salmon'
    @params['transcriptFasta'] = ''
    @params['transcriptFasta', 'description'] = 'give full path of transcript fasta file; in that case the build is ignored; if it comes from trinity assembly the gene-isoform associations will be extracted and used'
    @params['transcriptTypes'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA', 'long_noncoding', 'short_noncoding', 'pseudogene']
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
    @params['poly_x_min_len'] = '10'
    @params['poly_x_min_len','description'] = 'the minimum length to detect polyX in the read tail. 10 by default.'
    @params['length_required'] = '18'
    @params['length_required','description'] = 'reads shorter than length_required will be discarded.'
    @params['cmdOptionsFastp'] = ''
    ## additional commands
    @params['markDuplicates'] = true
    @params['markDuplicates', 'description'] = 'should duplicates be marked with sambamba. It is recommended for ChIP-seq and ATAC-seq data.'

    @params['mail'] = ""
    @modules = ["Aligner/Salmon", "QC/Flexbar", "QC/fastp", "Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    else
      @required_columns << 'fldMean' << 'fldSD'
    end
  end
  def set_default_parameters
    @params['paired'] = dataset_has_column?('Read2')
    if dataset_has_column?('strandMode')
      @params['strandMode'] = @dataset[0]['strandMode']
    end
  end
  def next_dataset
    {
      'Name'=>@dataset['Name'],
      'Count [File]'=>File.join(@result_dir, "#{@dataset['Name']}.txt"),
      'Species'=>@dataset['Species'],
      'refBuild'=>@params['refBuild'],
      'featureLevel'=>'isoform',
      'refFeatureFile'=>@params['refFeatureFile'],
      'strandMode'=>@params['strandMode'],
      'paired'=>@params['paired'],
      'Read Count'=>@dataset['Read Count'],
      'transcriptTypes'=>@params['transcriptTypes']
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppSalmon")
  end
end

if __FILE__ == $0
  run SalmonApp

end
