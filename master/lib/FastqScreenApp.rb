#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class FastqScreenApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'FastqScreen'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'QC'
@description =<<-EOS
Screen files for contaminations or ribosomal RNA content<br/>
<a target='_blank' href='http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/'>fastq_screen web site</a>
EOS
    @required_columns = ['Name','Read1','Read Count']
    @required_params = ['name']
    @params['cores'] = '8'
    @params['ram'] = '60'
    @params['scratch'] = '100'
    @params['name'] = 'FastqScreen_Result'
    @params['nReads'] = '100000'
    @params['paired'] = false
    @params['readFileToUse', 'description'] = 'Which reads should be used for analysis. Can be used to override "paired" parameter.'
    @params['readFileToUse'] = ['Read1', 'Read2', 'both']
    @params['nTopSpecies'] = '5'
    @params['minAlignmentScore'] = '-20'
    @params['minAlignmentScore', 'description'] = 'the alignment score for bowtie2; can be negative'
    #@params['confFile'] = {'select'=>''}
    #if defined?(FASTQSCREEN_CONF_DIR)
    #  Dir["#{FASTQSCREEN_CONF_DIR}/*.conf"].sort.select{|conf| File.file?(conf)}.each do |file|
    #    @params['confFile'][File.basename(file)] = File.basename(file)
    #  end
    #end

    @params['cmdOptions'] = "-k 10 --trim5 4 --trim3 4 --very-sensitive"
    
    # general
    # @params['trimAdapter', 'hr'] = true
    # @params['trimAdapter'] = true
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
    @params['mail'] = ""
    @modules = ["Tools/samtools", "Aligner/Bowtie2", "QC/fastp", "Tools/kraken", "QC/FastQScreen", "Dev/R", "Tools/Picard"]
    @inherit_columns = ["Order Id"]
  end
  def set_default_parameters
    @params['paired'] = dataset_has_column?('Read2')
  end
  def preprocess
    if @params['readFileToUse'] != 'Read1'
      @required_columns<<  'Read2'
    end
    if !dataset_has_column?('Read2')
      @params['readFileToUse'] = ['Read1']
    end
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
    run_RApp("EzAppFastqScreen")
  end
end


if __FILE__ == $0
  usecase = FastqscreenApp.new

  usecase.project = "p1001"
  usecase.user = "masa"

  # set user parameter
  # for GUI sushi
  #usecase.params['process_mode'].value = 'SAMPLE'
  #usecase.params['refBuild'] = 'TAIR10'
  #usecase.params['paired'] = true
  #usecase.params['cores'] = 2
  #usecase.params['node'] = 'fgcz-c-048'

  # also possible to load a parameterset csv file
  # mainly for CUI sushi
  usecase.parameterset_tsv_file = 'tophat_parameterset.tsv'
  #usecase.params['name'] = 'name'

  # set input dataset
  # mainly for CUI sushi
  usecase.dataset_tsv_file = 'tophat_dataset.tsv'

  # also possible to load a input dataset from Sushi DB
  #usecase.dataset_sushi_id = 1

  # run (submit to workflow_manager)
  usecase.run
  #usecase.test_run

end
