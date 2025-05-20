#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class PreqcApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Preqc'
    @analysis_category = 'QC'
    @description =<<-EOS
Preqc - Illumina read pre-assembly quality control and data exploration module within sga
<a href='https://github.com/jts/sga/wiki/Preqc'>https://github.com/jts/sga/wiki/Preqc</a>
EOS

    @required_columns = ['Name','Read1', 'Read2']
    @required_params = ['peMode', 'dustThreshold', 'example']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '100'
    @params['paired'] = true
    
        # preprocess options
    @params['peMode'] = ['1', '0', '2']
    @params['peMode', 'description'] = '0 - treat reads as single-end; 1 - reads are paired with the first read in R1.fastq and the second in R2.fastq; 2 - reads are paired and interleaved within a single file. Default is 1'
    @params['dustThreshold'] = '4.0'
    @params['dustThreshold', 'description'] = 'filter out reads that have a dust score higher than FLOAT (default: 4.0)'
    @params['example'] = ['human', 'bird', 'fish', 'oyster', 'snake', 'yeast']
    @params['example', 'description'] = 'pre-computed readset to be included in the final report'
        # trimming options
    # general
    #@params['trimAdapter', 'hr'] = true
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
    @params['mail'] = ""
    @modules = ["Dev/R", "QC/fastp", "Dev/Python/3.8.3"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'PreqcReport [Link]'=>File.join(@result_dir, "#{@dataset['Name']}.pdf"),
     'PreqcReport [File]'=>File.join(@result_dir, "#{@dataset['Name']}.pdf"),
     'PreqcOut [File]'=>File.join(@result_dir, "#{@dataset['Name']}.preqc"),
    }.merge(extract_columns(@inherit_tags))
  end
 def commands
    run_RApp("EzAppPreqc",conda_env: "gi_sga0.10.15")
  end 
end

if __FILE__ == $0
end
