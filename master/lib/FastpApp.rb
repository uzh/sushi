#!/usr/bin/env ruby
# encoding: utf-8
Version = '20200218-163350'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class FastpApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'fastp'
    @analysis_category = 'Prep'
    @description =<<-EOS
A tool designed to provide fast all-in-one preprocessing for FastQ files. This tool is developed in C++ with multithreading supported to afford high performance.
Refer to <a href='https://github.com/OpenGene/fastp'>https://github.com/OpenGene/fastp</a>
EOS
    @required_columns = ['Name','Read1']
    @required_params = ['paired']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '16'
    @params['scratch'] = '100'
    @params['paired'] = false
    @params['paired', 'description'] = 'either the reads are paired-ends or single-end'
    # trimming options
    # general
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
    @params['mail'] = ""
    @modules = ["QC/fastp"]
    @inherit_tags = ["Factor"]
  end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  end
  def set_default_parameters
    @params['paired'] = dataset_has_column?('Read2')
  end
  def next_dataset
   dataset =  {
      'Name'=>@dataset['Name'],
      'Read1 [File]' => File.join(@result_dir, "#{File.basename(@dataset['Read1'].to_s).gsub('fastq.gz','trimmed.fastq.gz')}"),
      'Adapters [File]' => File.join(@result_dir, "#{@dataset['Name']}_adapters.fa")
    }.merge(extract_columns(@inherit_tags))
    if @params['paired'] 
        dataset['Read2 [File]'] = File.join(@result_dir, "#{File.basename(@dataset['Read2'].to_s).gsub('fastq.gz','trimmed.fastq.gz')}")
    end
    dataset
  end
  def commands
    run_RApp("EzAppFastp")
  end
end

if __FILE__ == $0
  run FastpApp
end

