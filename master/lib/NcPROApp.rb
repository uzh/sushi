#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class NcPROApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'ncPRO_Report'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Count'
    @description =<<-EOS
Annotation and Profiling of ncRNAs in smallRNA-seq<br/>
Uses ncPRO-seq for a complete analysis of small-RNA-seq. ncPRO-seq considers performs quality assessment and quantiation.
It considers various species of small RNA.<br />
<a href='https://ncpro.curie.fr/index.html'>https://ncpro.curie.fr/index.html/</a>
ncPRO can only process single-end stranded RNA.
<strong>IMPORTANT: ncPRO does not run on the virtual nodes!!!!</strong>
EOS

    @required_columns = ['Name','Read1', 'Adapter1', 'Species']
    @required_params = ['name', 'cores', 'ram', 'scratch']

    @params['cores'] = '8'
    @params['ram'] = '40'
    @params['scratch'] = '100'
    @params['refBuild'] = ['Mus_musculus/UCSC/mm10', 'Homo_sapiens/UCSC/hg19', 'Rattus_norvegicus/UCSC/rn5']
    @params['paired'] = false
    @params['paired', 'description'] = 'whether the reads are paired end; must be false since ncPRO-seq does not support paired-end'
    @params['name'] = 'ncPRO_Result'
    
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
    @modules = ["Aligner/Bowtie", "QC/Flexbar", "Tools/ncPROseq", "QC/Trimmomatic", "Dev/R"]
  end
  def preprocess
    if @params['paired']
      @required_columns<<  'Read2'
    end
  end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, 'ncpro/report.html')
    {'Name'=>@params['name'],
     'Species'=>(dataset = @dataset.first and dataset['Species']),
     'refBuild'=>@params['refBuild'],
     'Report [File]'=>report_file,
     'Html [Link]'=>report_link,
     'TrimCounts [Link]'=>File.join(report_file, 'trimCounts-barplot.png')
    }
  end
  def commands
    run_RApp("EzAppNcpro")
  end
end

if __FILE__ == $0
  run NcPROApp

end
