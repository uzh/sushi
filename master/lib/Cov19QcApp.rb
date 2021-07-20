#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class Cov19QcApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'Cov19QC'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'QC'
    @description =<<-EOS 
Simple quality control for Cov19 data<br/>
EOS
    @required_columns = ['Name','Read1']
    @required_params = ['name', 'paired']
    @params['cores'] = [8, 1, 2, 4, 8]
    @params['ram'] = [15, 30]
    @params['ram', 'description'] = "GB"
    @params['scratch'] = [100, 10, 50, 100, 200]
    @params['scratch', 'description'] = "GB"
    @params['paired'] = false
    @params['refBuild'] = "Human_coronavirus_2019/ncbi/MN908947.V3"
    @params['name'] = "CovidQC"
    @params['readsUsed']="500000"
    @params['Adapter1']="CTGTCTCTTATACACATCT"
    @params['cmdOptions'] = ""

    # trimming options
    # general
    @params['trimAdapter'] = true
    # Fastp
    ## trimming
    @params['trim_front1'] = '4'
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
    @params['average_qual'] = '10'
    @params['average_qual','description'] = 'if one read average quality score <avg_qual>, then this read/pair is discarded. Default 0 means no requirement'
    @params['max_len1'] = '0'
    @params['max_len1','description'] = 'if read1 is longer than max_len1, then trim read1 at its tail to make it as long as max_len1. Default 0 means no limitation.'
    @params['max_len2'] = '0'
    @params['max_len2','description'] = 'if read1 is longer than max_len2, then trim read2 at its tail to make it as long as max_len1. Default 0 means no limitation.'
    @params['poly_x_min_len'] = '10'
    @params['poly_x_min_len','description'] = 'the minimum length to detect polyX in the read tail. 10 by default.'
    @params['length_required'] = '30'
    @params['length_required','description'] = 'reads shorter than length_required will be discarded.'
    
    @params['mail'] = ""
    @modules = ["Aligner/Bowtie2", "Dev/R", "QC/fastp", "Tools/samtools", "Tools/bcftools", "Dev/jdk", "Tools/Picard"]
    @inherit_tags = ["B-Fabric"]
  end
 def set_default_parameters
    @params['paired'] = dataset_has_column?('Read2')
  end
  def preprocess
    if @params['paired']
      @required_columns<<  'Read2'
    end
  end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@params['name'],
     'Report [File]'=>report_file,
     'Html [Link]'=>report_link
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppCov19QC")
  end
end
