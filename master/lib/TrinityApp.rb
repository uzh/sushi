#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class TrinityApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Trinity'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Assemble'
    @description =<<-EOS
    RNA-seq de novo Assembly<br/>
<a href='https://github.com/trinityrnaseq/trinityrnaseq/wiki'>Trinity</a><br/>
    EOS
    @required_columns = ['Name', 'Read1', 'Species']
    @required_params = ['name']
    # optional params
    @params['cores'] = '12'
    @params['ram'] = '220'
    @params['scratch'] = '1000'
    @params['name'] = "Trinity_Assembly"
    @params['paired'] = false
    @params['strandMode'] = ['both', 'sense', 'antisense']
    
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
    @params['markDuplicates'] = true
    @params['markDuplicates', 'description'] = 'should duplicates be marked with picard. It is recommended for ChIP-seq and ATAC-seq data.'

    @params['mail'] = ""
    @modules = ["Assembly/Trinity", "QC/fastp", "QC/Flexbar", "Aligner/Bowtie2", "Aligner/Salmon", "Dev/R"]
  end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  end
  def set_default_parameters
    @params['paired'] = dataset_has_column?('Read2')
    if dataset_has_column?('strandMode')
      @params['strandMode'] = @dataset[0]['strandMode']
    end
  end
  def next_dataset
    {'Name'=>@params['name'],
     'Fasta [File]'=>File.join(@result_dir, "#{@params['name']}.fasta"),
	  'Stats [File]'=>File.join(@result_dir, "assembly_stats.txt"),
	  'Abundance Counts [File]'=>File.join(@result_dir, "abundance_counts.txt"),
	  'Abundance TPM [File]'=>File.join(@result_dir, "abundance_TPM.txt"),
	  'Abundance TMM [File]'=>File.join(@result_dir, "abundance_TMM.txt"),
	  'ExN50 [File]'=>File.join(@result_dir, "ExN50_stats.txt")
    }
  end
  def commands
    run_RApp('EzAppTrinity')
  end
end

if __FILE__ == $0
  #usecase = EdgeRApp.new

  usecase.project = "p1001"
  usecase.user = 'masamasa'

  # set user parameter
  # for GUI sushi
  #usecase.params['process_mode'].value = 'SAMPLE'
  usecase.params['refBuild'] = 'mm10'
  usecase.params['paired'] = true
  usecase.params['strandMode'] = 'both'
  usecase.params['cores'] = 8
  usecase.params['node'] = 'fgcz-c-048'

  # also possible to load a parameterset csv file
  # mainly for CUI sushi
  #usecase.parameterset_tsv_file = 'tophat_parameterset.tsv'
  #usecase.parameterset_tsv_file = 'test.tsv'

  # set input dataset
  # mainly for CUI sushi
  #usecase.dataset_tsv_file = 'tophat_dataset.tsv'

  # also possible to load a input dataset from Sushi DB
  usecase.dataset_sushi_id = 3

  # run (submit to workflow_manager)
  usecase.run
  #usecase.test_run

end
