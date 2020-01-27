#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class BWAApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'BWA'
    @analysis_category = 'Map'
    @description =<<-EOS
    Burrows-Wheeler Aligner<br/>
<a href='http://bio-bwa.sourceforge.net/'>BWA</a><br/>
    EOS
    @required_columns = ['Name','Read1','Species']
    @required_params = ['refBuild','paired']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '16'
    @params['scratch'] = '100'
    @params['refBuild'] = ref_selector
    @params['refBuild', 'description'] = 'the genome refBuild and annotation to use as reference. If human variant calling is the main goal, please use hg_19_karyotypic.'
    @params['paired'] = false
    @params['algorithm'] = ['mem', 'aln', 'bwasw']
    @params['cmdOptions'] = ''
    
    # trimming options
    # general
    @params['trimAdapter'] = true
    # Fastp
    ## trimming
    @param['reads_to_process'] = '0'
    @param['reads_to_process','description'] = 'specify how many reads/pairs to be processed. Default 0 means process all reads.'
    @param['trim_front'] = '0'
    @param['trim_front','description'] = 'trimming how many bases in front for read1 (and read2), default is 0.'
    @param['trim_tail'] = '0'
    @param['trim_tail','description'] = 'trimming how many bases in tail for read1 (and read2), default is 0.'
    @param['cut_front'] = '0'
    @param['cut_front','description'] = 'move a sliding window from front (5p) to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.'
    @param['cut_tail'] = '0'
    @param['cut_tail','description'] = 'move a sliding window from tail (3p) to front, drop the bases in the window if its mean quality < threshold, stop otherwise.'
    @param['cut_right'] = '0'
    @param['cut_right','description'] = 'move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, and then stop.'
    @param['average_qual'] = '0'
    @param['average_qual','description'] = 'if one read average quality score <avg_qual>, then this read/pair is discarded. Default 0 means no requirement'
    @param['max_len1'] = '0'
    @param['max_len1','description'] = 'if read1 is longer than max_len1, then trim read1 at its tail to make it as long as max_len1. Default 0 means no limitation. If two reads are present, the same will apply to read2.'
    @param['poly_x_min_len'] = '10'
    @param['poly_x_min_len','description'] = 'the minimum length to detect polyX in the read tail. 10 by default.'
    @param['length_required'] = '15'
    @param['length_required','description'] = 'reads shorter than length_required will be discarded, default is 15.'
    @param['length_limit'] = '0'
    @param['length_limit','description'] = 'reads longer than length_limit will be discarded, default 0 means no limitation.'
    ## additional commands
    @params['cmdOptionsFastp'] = ''
    
    @params['mail'] = ""
    @modules = ["Tools/samtools", "Aligner/BWA", "QC/Flexbar", "QC/Trimmomatic", "Dev/R", "Tools/sambamba"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
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
    {'Name'=>@dataset['Name'],
     'BAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam"),
     'BAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam.bai"),
        'IGV Starter [Link]'=>File.join(@result_dir, "#{@dataset['Name']}-igv.jnlp"),
        'Species'=>@dataset['Species'],
        'refBuild'=>@params['refBuild'],
        'paired'=>@params['paired'],
        'refFeatureFile'=>@params['refFeatureFile'],
        'strandMode'=>@params['strandMode'],
        'Read Count'=>@dataset['Read Count'],
        'IGV Starter [File]'=>File.join(@result_dir, "#{@dataset['Name']}-igv.jnlp"),
        'IGV Session [File]'=>File.join(@result_dir, "#{@dataset['Name']}-igv.xml"),
        'PreprocessingLog [File]'=>File.join(@result_dir, "#{@dataset['Name']}_preprocessing.log")
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppBWA")
  end
end

if __FILE__ == $0

end
