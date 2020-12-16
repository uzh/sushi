#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-093702'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class Bowtie2TranscriptomeApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Bowtie2TranscriptomeApp'
    @analysis_category = 'Map'
    @description =<<-EOS
Aligns the reads to a transcriptome only index. The index is built from the gtf file specified. For alignments to
genome or joined see the other mappers.
Mapping to transcriptome is useful, e.g. for Ribo-Seq data.
For expression estimation of RNA-seq data we recommend using the RSEM counting App or STAR+featureCounts.
For bowtie2 options see:<br/>
<a href='http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml'>http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml/</a>
EOS

    @required_columns = ['Name','Read1','Species']
    @required_params = ['refBuild','paired']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '15'
    @params['scratch'] = '100'
    @params['refBuild'] = ref_selector
    @params['refBuild', 'description'] = 'the genome refBuild and annotation to use as reference. If human variant calling is the main goal, please use hg_19_karyotypic.'
    @params['paired'] = false
    @params['paired', 'description'] = 'whether the reads are paired end; if false then only Read1 is considered even if Read2 is available.'
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['refFeatureFile'] = 'genes.gtf'
    @params['cmdOptions'] = '--no-unal'
    @params['cmdOptions', 'description'] = 'specify the commandline options for bowtie2; do not specify any option that is already covered by the dedicated input fields'
    
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
    @params['max_len2'] = '0'
    @params['max_len2','description'] = 'if read1 is longer than max_len2, then trim read2 at its tail to make it as long as max_len1. Default 0 means no limitation.'
    @params['poly_x_min_len'] = '10'
    @params['poly_x_min_len','description'] = 'the minimum length to detect polyX in the read tail. 10 by default.'
    @params['length_required'] = '18'
    @params['length_required','description'] = 'reads shorter than length_required will be discarded.'
    ## additional commands
    @params['cmdOptionsFastp'] = ''
    
    @params['mail'] = ""
    @modules = ["Tools/samtools", "Aligner/Bowtie2", "QC/fastp", "Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
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
    {'Name'=>@dataset['Name'],
     'trBAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam"),
     'trBAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam.bai"),
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'paired'=>@params['paired'],
     'refFeatureFile'=>@params['refFeatureFile'],
     'strandMode'=>@params['strandMode'],
     'Read Count'=>@dataset['Read Count']
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppBowtie2Transcriptome")
  end
end
