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
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '100'
    @params['refBuild'] = ref_selector
    @params['refBuild', 'description'] = 'the genome refBuild and annotation to use as reference.'
    @params['paired'] = false
    @params['paired', 'description'] = 'whether the reads are paired end; if false then only Read1 is considered even if Read2 is available.'
    @params['deduplicate'] = true
    @params['deduplicate', 'description'] = 'Removal of PCR duplicates. Not recommended for RRBS-Seq (Enzymatic digestion)'
    @params['cmdOptions'] = '--bowtie2'
    @params['cmdOptions', 'description'] = 'specify the commandline options for bismark; do not specify any option that is already covered by the dedicated input fields'
    @params['trimAdapter'] = true
    @params['trimAdapter', 'description'] = 'if adapters should be trimmed'
    @params['trimLeft'] = 0
    @params['trimLeft', 'description'] = 'fixed trimming at the "left" i.e. 5-prime end of the read'
    @params['trimRight'] = 0
    @params['trimRight', 'description'] = 'fixed trimming at the "right" i.e. 3-prime end of the read'
    @params['minTailQuality'] = 0
    @params['minTailQuality', 'description'] = 'if above zero, then reads are trimmed as soon as 4 consecutive bases have lower mean quality'
    @params['minAvgQuality'] = 10
    @params['minReadLength'] = 20
    @params['specialOptions'] = ''
    @params['specialOptions', 'description'] = 'special unsupported options that the R wrapper may support, format: <key>=<value>'
    @params['mail'] = ""
    @modules = ["Tools/samtools", "Aligner/Bowtie2", "Aligner/Bismark", "QC/Flexbar", "QC/Trimmomatic", "Dev/R", "Tools/sambamba"]
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
