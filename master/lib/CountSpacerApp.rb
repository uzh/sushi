#!/usr/bin/env ruby
# encoding: utf-8
Version = '20200506-061900'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class CountSpacerApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'CountSpacer'
    @analysis_category = 'Misc'
    @description =<<-EOS
QC Tool for sgRNA libraries.
<a href='https://github.com/fengzhanglab/Screening_Protocols_manuscript'>https://github.com/fengzhanglab/Screening_Protocols_manuscript</a>
EOS

    @required_columns = ['Name','Read1']
    @required_params = ['dictPath']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '40'
    @params['scratch'] = '100'

    @params['name'] = 'CountSpacer_Result'
    @params['dictPath'] = ''
    @params['dictPath'] = {'select'=>''}
    Dir["/srv/GT/databases/GEML/sgRNA_Libs/*"].sort.select{|lib| File.directory?(lib)}.each do |dir|
      @params['dictPath'][File.basename(dir)] = File.basename(dir)
    end
    @params['leftPattern'] = ''
    @params['leftPattern', 'description'] = 'short patterns < 8bp could could cause misleading results'
    @params['rightPattern'] = ''
    @params['rightPattern', 'description'] = 'short patterns < 8bp could could cause misleading results'
    @params['maxMismatch'] = 1
    @params['maxMismatch', 'description'] = 'number of allowed mismatches for pattern search'

    @params['specialOptions'] = ''
    @params['specialOptions', 'description'] = 'special unsupported options that the R wrapper may support, format: <key>=<value>'
    
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
    @params['average_qual'] = '10'
    @params['average_qual','description'] = 'if one read average quality score <avg_qual>, then this read/pair is discarded. Default 0 means no requirement'
    @params['max_len1'] = '0'
    @params['max_len1','description'] = 'if read1 is longer than max_len1, then trim read1 at its tail to make it as long as max_len1. Default 0 means no limitation. If two reads are present, the same will apply to read2.'
    @params['max_len2'] = '0'
    @params['max_len2','description'] = 'if read1 is longer than max_len2, then trim read2 at its tail to make it as long as max_len1. Default 0 means no limitation.'
    @params['poly_x_min_len'] = '10'
    @params['poly_x_min_len','description'] = 'the minimum length to detect polyX in the read tail. 10 by default.'
    @params['length_required'] = '19'
    @params['length_required','description'] = 'reads shorter than length_required will be discarded.'
    @params['cmdOptionsFastp'] = ''
    ## additional commands
    @params['mail'] = ""
    @modules = ["QC/fastp", "Dev/R", "Aligner/Bowtie"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def preprocess
    #if @params['paired']
    #  @required_columns << 'Read2'
    #end
  end
 def set_default_parameters
    #@params['paired'] = dataset_has_column?('Read2')
  end
  def next_dataset
    report_file = File.join(@result_dir,"#{@dataset['Name']}")
    report_link = File.join(report_file, '00index.html')
    count_file = File.join(report_file, "#{@dataset['Name']}-sgRNA_counts.xlsx")
    {'Name'=>@dataset['Name'],
     'Report [File]'=>report_file,
     'Html [Link]'=>report_link,
     'Count [File]'=>count_file, 
     'Species'=>@dataset['Species'],
     #'refBuild'=>'',
     #'featureLevel'=>'smRNA',
     #'refFeatureFile'=>'',
     'Read Count'=>@dataset['Read Count'] 
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppCountSpacer")
  end
end

if __FILE__ == $0
  run CountSpacerApp
  #usecase = Bowtie2App.new

  #usecase.project = "p1001"
  #usecase.user = 'masamasa'

  # set user parameter
  # for GUI sushi
  #usecase.params['process_mode'].value = 'SAMPLE'
  #usecase.params['refBuild'] = 'mm10'
  #usecase.params['paired'] = true
  #usecase.params['strandMode'] = 'both'
  #usecase.params['cores'] = 8
  #usecase.params['node'] = 'fgcz-c-048'

  # also possible to load a parameterset csv file
  # mainly for CUI sushi
  #usecase.parameterset_tsv_file = 'tophat_parameterset.tsv'
  #usecase.parameterset_tsv_file = 'test.tsv'

  # set input dataset
  # mainly for CUI sushi
  #usecase.dataset_tsv_file = 'tophat_dataset.tsv'

  # also possible to load a input dataset from Sushi DB
  #usecase.dataset_sushi_id = 3

  # run (submit to workflow_manager)
  #usecase.run
  #usecase.test_run

end
