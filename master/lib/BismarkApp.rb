#!/usr/bin/env ruby
# encoding: utf-8
Version = '20180408-135252'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class BismarkApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Bismark'
    @analysis_category = 'Map'
    @description =<<-EOS
A tool to map bisulfite converted sequence reads and determine cytosine methylation states<br/>
<a href='http://www.bioinformatics.babraham.ac.uk/projects/bismark/'>http://www.bioinformatics.babraham.ac.uk/projects/bismark</a>
EOS

    @required_columns = ['Name','Read1','Species']
    @required_params = ['refBuild','paired']
    # optional params
    @params['cores'] = ['8', '12', '16']
    @params['cores', "context"] = "slurm"
    @params['ram'] = '30'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '200'
    @params['scratch', "context"] = "slurm"
    @params['refBuild'] = ref_selector
    @params['refBuild', 'description'] = 'the genome refBuild and annotation to use as reference.'
    @params['refBuild', "context"] = "referfence genome assembly"
    @params['paired'] = false
    @params['paired', 'description'] = 'whether the reads are paired end; if false then only Read1 is considered even if Read2 is available.'
    @params['paired', "context"] = "Bismark"
    @params['deduplicate'] = true
    @params['deduplicate', 'description'] = 'Removal of PCR duplicates. Not recommended for RRBS-Seq (Enzymatic digestion)'
    @params['deduplicate', "context"] = "Bismark"
    @params['cmdOptions'] = '--bowtie2'
    @params['cmdOptions', 'description'] = 'specify the commandline options for bismark; do not specify any option that is already covered by the dedicated input fields'
    @params['cmdOptions', "context"] = "Bismark"
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
    @params['generateBigWig'] = false
    @params['generateBigWig', "context"] = "Bismark"
    @params['EM_QC'] = false
    @params['EM_QC','description'] = 'Generate Boxplot for known methylated and unmethylated sites of Lambda/pUC19 control'
    @params['EM_QC', "context"] = "Bismark"
    @params['mail'] = ""
    @modules = ["Tools/samtools", "Aligner/Bowtie2", "Aligner/Bismark", "QC/fastp", "Dev/R", "Dev/Perl"]
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
   dataset = {
     'Name'=>@dataset['Name'],
     'BAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam"),
     'BAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam.bai"),
     'TxtReport [File]'=>File.join(@result_dir, "#{@dataset['Name']}.report.txt"),
     'M-Bias_R1 [File]'=>File.join(@result_dir, "#{@dataset['Name']}.M-bias_R1.png"),
     'M-Bias_R2 [File]'=>File.join(@result_dir, "#{@dataset['Name']}.M-bias_R2.png"),
     'CpG_Context [File]'=>File.join(@result_dir, "#{@dataset['Name']}.CpG_context.txt"),
     'COV [File]'=>File.join(@result_dir, "#{@dataset['Name']}.gz.bismark.cov.gz"),
     'BedGraph [File]'=>File.join(@result_dir, "#{@dataset['Name']}.gz"),
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'paired'=>@params['paired'],
     'Read Count'=>@dataset['Read Count'],
     'PreprocessingLog [File]'=>File.join(@result_dir, "#{@dataset['Name']}_preprocessing.log")
    }.merge(extract_columns(@inherit_tags))
    if @params['generateBigWig']
       dataset['BigWig [File]'] = File.join(@result_dir, "#{@dataset['Name']}_Cov.bw")
    end
   if @params['EM_QC']
       dataset['EM_QC [Link]'] = File.join(@result_dir, "#{@dataset['Name']}_Meth.png")
       dataset['EM_QC [File]'] = File.join(@result_dir, "#{@dataset['Name']}_Meth.png")
    end
    dataset
  end
  def commands
    run_RApp("EzAppBismark")
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
