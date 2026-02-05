#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class FastqScreen10xApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'FastqScreen10x'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'QC'
@description =<<-EOS
Screen files for contaminations or ribosomal RNA content<br/>
<a target='_blank' href='http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/'>fastq_screen web site</a>
EOS
    @required_columns = ['Name','RawDataDir', 'Read Count']
    @required_params = ['name', 'paired']
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '80'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '300'
    @params['scratch', "context"] = "slurm"
    @params['paired'] = false
    @params['paired', "context"] = "FastqScreen10x"
    @params['readFileToUse'] = ['Read2']
    @params['readFileToUse', "context"] = "FastqScreen10x"
    @params['name'] = 'FastqScreen_Result'
    @params['nReads'] = '100000'
    @params['nReads', "context"] = "FastqScreen10x"
    @params['nTopSpecies'] = '5'
    @params['nTopSpecies', "context"] = "FastqScreen10x"
    @params['minAlignmentScore'] = '-20'
    @params['minAlignmentScore', 'description'] = 'the alignment score for bowtie2; can be negative'
    @params['minAlignmentScore', "context"] = "FastqScreen10x"
    @params['cmdOptions'] = "-k 10 --trim5 4 --trim3 4 --very-sensitive"
    @params['cmdOptions', "context"] = "FastqScreen10x"
    
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
    @params['mail'] = ""
    
    @modules = ["Aligner/BWA", "Tools/samtools", "Aligner/Bowtie2", "QC/fastp", "Tools/kraken", "QC/FastQScreen", "Tools/Picard"]
    @inherit_columns = ["Order Id"]
  end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@params['name'],
     'Report [File]'=>report_file,
     'Html [Link]'=>report_link,
    }.merge(extract_columns(colnames: @inherit_columns))
  end
  def commands
    commands = "export LANG=en_US.UTF-8 \n"
    commands << run_RApp("EzAppFastqScreen_10x")
  end
end
