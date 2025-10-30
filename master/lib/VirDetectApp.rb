#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-100301'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class VirDetectApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'VirDetect'
    @analysis_category = 'Metagenomics'
    @description =<<-EOS
Virome analysis for the diagnostic use in veterinary virology<br/>
<a href='http://www.vetvir.uzh.ch/en/Research/Virology/Molecular-and-Clinical-Veterinary-Virology/Virome-Analysis.html'>Virome Analysis Group</a>, The Institute of Virology, UZH<br>
<a href='https://www.ncbi.nlm.nih.gov/pubmed/25009045'>Reference paper</a>, base on which this pileline is implemented.
EOS

    @required_columns = ['Name','Read1','Species']
    @required_params = ['virBuild','hostBuild','paired']
    @params['virBuild'] = ref_selector
    @params['virBuild', 'description'] = 'the viral reference database to use.'
    @params['virBuild', "context"] = "VirDetect"
    @params['hostBuild'] = ref_selector
    @params['hostBuild', 'description'] = 'the none-human host genome to use as reference for removal of contamination'
    @params['hostBuild', "context"] = "VirDetect"
    @params['paired'] = false
    @params['paired', 'description'] = 'whether the reads are paired end; if false then only Read1 is considered even if Read2 is available.'
    @params['paired', "context"] = "VirDetect"
    # optional params
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '60'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '150'
    @params['scratch', "context"] = "slurm"
    @params['cmdOptionsHost'] = '--very-sensitive'
    @params['cmdOptionsHost', 'description'] = 'specify the commandline options for bowtie2 during mapping to contaminated genome(s); do not specify any option that is already covered by the dedicated input fields'
    @params['cmdOptionsHost', "context"] = "VirDetect"
    @params['cmdOptions'] = '-a --very-sensitive --no-mixed --no-discordant -X 1000'
    @params['cmdOptions', 'description'] = 'specify the commandline options for bowtie2 during mapping to the viral reference database; do not specify any option that is already covered by the dedicated input fields'
    @params['cmdOptions', "context"] = "VirDetect"
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
    
    @modules = ["Tools/samtools", "Aligner/Bowtie2", "QC/fastp", "Tools/BEDTools", "Dev/R"]
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
     'Species'=>@dataset['Species'],
     'virBuild'=>@params['virBuild'],
     'hostBuild'=>@params['hostBuild'],
     'paired'=>@params['paired'],
     'Read Count'=>@dataset['Read Count'],
     'OutDir [File]'=>File.join(@result_dir, "#{@dataset['Name']}"),
     'OutReport [Link]'=>File.join(@result_dir, "#{@dataset['Name']}", "#{@dataset['Name']}.html")
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppVirDetect")
  end
end

if __FILE__ == $0
  run Bowtie2App
  #usecase = VirDetectApp.new

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
