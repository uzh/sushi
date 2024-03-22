#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class UnicyclerApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Unicycler'
    @analysis_category = 'Assemble'
    @description =<<-EOS
Unicycler microbial genome assembler
<a href='https://github.com/rrwick/Unicycler'>https://github.com/rrwick/Unicycler</a>
EOS

    @required_columns = ['Name','Read1']
    @required_params = ['cores', 'ram', 'scratch']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '100'
    @params['paired'] = true
    @params['long'] = ["NO", "YES"]
    @params['long', 'description'] = "Do you have long reads with which you would like to do a hybrid assembly"
    @params['pathToLong'] = ''
    @params['pathToLong', 'description'] = "Specify the path to long reads, make sure the file has the right permissions."
    @params['mode'] = ['normal', 'conservative', 'bold']
    @params['mode', 'description'] = "Unicycler can be run in three modes. Default is normal. See documentation for more clarification."
    @params['unicyclerOpt'] = ''
    @params['unicyclerOpt', 'description'] = 'Predefine any options not already specified by default. By default is empty for genome assembly'
    
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
    #@params['markDuplicates'] = true
    #@params['markDuplicates', 'description'] = 'should duplicates be marked with picard. It is recommended for ChIP-seq and ATAC-seq data.'

    @params['mail'] = ""
    @modules = ["Assembly/SPAdes", "QC/fastp", "Dev/R", "Assembly/Unicycler"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'Assembly [File]'=>File.join(@result_dir, "assembly.fasta"),
     'Graph [File]'=>File.join(@result_dir, "assembly.gfa"),
     'Log [File]'=>File.join(@result_dir, "unicycler.log"),
     'Species'=>@dataset['Species'],
     'Read Count'=>@dataset['Read Count'],
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppUnicycler")
  end
end

if __FILE__ == $0
end
