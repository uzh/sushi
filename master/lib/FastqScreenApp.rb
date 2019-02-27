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
    @required_columns = ['Name','Read1']
    @required_params = ['name', 'paired']
    @params['cores'] = '8'
    @params['ram'] = '40'
    @params['scratch'] = '100'
    @params['paired'] = false
    @params['name'] = 'FastqScreen_Result'
    #@params['trimAdapter'] = false
    @params['trimLeft'] = 0
    @params['trimRight'] = 0
    @params['minTailQuality'] = 10
    @params['minAvgQuality'] = 20
    @params['minReadLength'] = 30
    @params['nReads'] = '100000'
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
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Aligner/BWA", "Tools/samtools", "Aligner/Bowtie2", "QC/Trimmomatic", "QC/Flexbar", "QC/FastQScreen", "Dev/R", "Tools/sambamba"]
  end
 def set_default_parameters
    @params['paired'] = dataset_has_column?('Read2')
  end
  def preprocess
    if @params['paired']
      @required_columns<<  'Read2'
    end
  end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@params['name'],
     'Report [File]'=>report_file,
     'Html [Link]'=>report_link,
    }
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
