#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-093644'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class Bowtie2App < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Bowtie2'
    @analysis_category = 'Map'
    @description =<<-EOS
Fast and sensitive read alignment. Supports local and end-to-end mode<br/>
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
    @params['max_len1','description'] = 'if read1 is longer than max_len1, then trim read1 at its tail to make it as long as max_len1. Default 0 means no limitation.'
    @params['max_len2'] = '0'
    @params['max_len2','description'] = 'if read1 is longer than max_len2, then trim read2 at its tail to make it as long as max_len1. Default 0 means no limitation.'
    @params['poly_x_min_len'] = '10'
    @params['poly_x_min_len','description'] = 'the minimum length to detect polyX in the read tail. 10 by default.'
    @params['length_required'] = '18'
    @params['length_required','description'] = 'reads shorter than length_required will be discarded.'
    ## additional commands
    @params['cmdOptionsFastp'] = ''
    @params['markDuplicates'] = true
    @params['markDuplicates', 'description'] = 'should duplicates be marked with sambamba. It is recommended for ChIP-seq and ATAC-seq data.'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Tools/samtools", "Aligner/Bowtie2", "QC/Flexbar", "QC/fastp", "Dev/R", "Tools/sambamba", "Tools/Picard"]
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
     'IGV [Link]'=>File.join(@result_dir, "#{@dataset['Name']}-igv.html"),
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'paired'=>@params['paired'],
     'Read Count'=>@dataset['Read Count'],
     'IGV [File]'=>File.join(@result_dir, "#{@dataset['Name']}-igv.html"),
     'PreprocessingLog [File]'=>File.join(@result_dir, "#{@dataset['Name']}_preprocessing.log"),
     'Bowtie2Log [File]'=>File.join(@result_dir, "#{@dataset['Name']}_bowtie2.log")
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppBowtie2")
  end
end

if __FILE__ == $0
  run Bowtie2App
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
