#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class KrakenApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Kraken'
    @analysis_category = 'Metagenomics'
    @description =<<-EOS
Kraken taxonomic sequence classification system
<a href='https://ccb.jhu.edu/software/kraken2/index.shtml'>https://ccb.jhu.edu/software/kraken2/index.shtml</a>
EOS

    @required_columns = ['Name','Read1']
    @required_params = ['paired', 'krakenDBOpt']
    # optional params
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '100'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '500'
    @params['scratch', "context"] = "slurm"
    @params['paired'] = false
    @params['paired', "context"] = "Kraken"
    
    @params['krakenDBOpt'] = ['Standard', 'Viral', 'PlusPF', 'core_nt','k2_standard']
    @params['krakenDBOpt', 'description'] = 'kraken database options (https://benlangmead.github.io/aws-indexes/k2). Default is Standard'
    @params['krakenDBOpt', "context"] = "Kraken"
    @params['krakenConfidenceOpt'] = '0.0'
    @params['krakenConfidenceOpt', 'description'] = 'Confidence score threshold, between 0 and 1'
    @params['krakenConfidenceOpt', "context"] = "Kraken"
    @params['krakenPhredOpt'] = '0'
    @params['krakenPhredOpt', 'description'] = 'Phred score threshold'
    @params['krakenPhredOpt', "context"] = "Kraken"
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'specify other commandline options for kraken; do not specify any option that is already covered by the dedicated input fields'
    @params['cmdOptions', "context"] = "Kraken"
    
        # trimming options
    # general
    @params['trimAdapter', 'hr'] = true
    @params['trimAdapter'] = true
    @params['trimAdapter', "context"] = "OpenGene/fastp"
    # Fastp
    ## trimming
    @params['trim_front1'] = '0'
    @params['trim_front1','description'] = 'trimming how many bases in front for read1 (and read2), default is 0.'
    @params['trim_front1', "context"] = "OpenGene/fastp"
    @params['trim_tail1'] = '0'
    @params['trim_tail1','description'] = 'trimming how many bases in tail for read1 (and read2), default is 0.'
    @params['trim_tail1', "context"] = "OpenGene/fastp"
    @params['cut_front'] = false
    @params['cut_front','description'] = 'move a sliding window from front (5p) to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.'
    @params['cut_front', "context"] = "OpenGene/fastp"
    @params['cut_front_window_size'] = '4'
    @params['cut_front_window_size', "context"] = "OpenGene/fastp"
    @params['cut_front_mean_quality'] = '20'
    @params['cut_front_mean_quality', "context"] = "OpenGene/fastp"
    @params['cut_tail'] = false
    @params['cut_tail','description'] = 'move a sliding window from tail (3p) to front, drop the bases in the window if its mean quality < threshold, stop otherwise.'
    @params['cut_tail', "context"] = "OpenGene/fastp"
    @params['cut_tail_window_size'] = '4'
    @params['cut_tail_window_size', "context"] = "OpenGene/fastp"
    @params['cut_tail_mean_quality'] = '20'
    @params['cut_tail_mean_quality', "context"] = "OpenGene/fastp"
    @params['cut_right'] = false
    @params['cut_right','description'] = 'move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, and then stop.'
    @params['cut_right', "context"] = "OpenGene/fastp"
    @params['cut_right_window_size'] = '4'
    @params['cut_right_window_size', "context"] = "OpenGene/fastp"
    @params['cut_right_mean_quality'] = '20'
    @params['cut_right_mean_quality', "context"] = "OpenGene/fastp"
    @params['average_qual'] = '0'
    @params['average_qual','description'] = 'if one read average quality score <avg_qual>, then this read/pair is discarded. Default 0 means no requirement'
    @params['average_qual', "context"] = "OpenGene/fastp"
    @params['max_len1'] = '0'
    @params['max_len1','description'] = 'if read1 is longer than max_len1, then trim read1 at its tail to make it as long as max_len1. Default 0 means no limitation. If two reads are present, the same will apply to read2.'
    @params['max_len1', "context"] = "OpenGene/fastp"
    @params['max_len2'] = '0'
    @params['max_len2','description'] = 'if read1 is longer than max_len2, then trim read2 at its tail to make it as long as max_len1. Default 0 means no limitation.'
    @params['max_len2', "context"] = "OpenGene/fastp"
    @params['poly_x_min_len'] = '10'
    @params['poly_x_min_len','description'] = 'the minimum length to detect polyX in the read tail. 10 by default.'
    @params['poly_x_min_len', "context"] = "OpenGene/fastp"
    @params['length_required'] = '18'
    @params['length_required','description'] = 'reads shorter than length_required will be discarded.'
    @params['length_required', "context"] = "OpenGene/fastp"
    @params['cmdOptionsFastp'] = ''
    @params['cmdOptionsFastp', "context"] = "OpenGene/fastp"
    ## additional commands
    @params['markDuplicates'] = true
    @params['markDuplicates', 'description'] = 'should duplicates be marked with picard. It is recommended for ChIP-seq and ATAC-seq data.'
    @params['markDuplicates', "context"] = "Picard"
    
    @params['mail'] = ""
    @modules = ["QC/fastp", "Dev/R", "Tools/kraken", "Tools/Krona"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'KronaReport [Link]'=>File.join(@result_dir, "#{@dataset['Name']}.html"),
     'KrakenOut [File]'=>File.join(@result_dir, "#{@dataset['Name']}.txt"),
     'KrakenReport [File]'=>File.join(@result_dir, "#{@dataset['Name']}.report.txt"),
     'KronaOutDir [File]'=>File.join(@result_dir, "#{@dataset['Name']}.html.files"),
     'KronaOut [File]'=>File.join(@result_dir, "#{@dataset['Name']}.html"),
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppKraken")
  end
end

if __FILE__ == $0
end
