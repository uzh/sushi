#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-094409'

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
    @params['trimAdapter'] = true
    @params['trimAdapter', 'description'] = 'if adapters should be trimmed'
    @params['trimLeft'] = 0
    @params['trimLeft', 'description'] = 'fixed trimming at the "left" i.e. 5-prime end of the read'
    @params['trimRight'] = 0
    @params['trimRight', 'description'] = 'fixed trimming at the "right" i.e. 3-prime end of the read'
    @params['minTailQuality'] = 0
    @params['minTailQuality', 'description'] = 'if above zero, then reads are trimmed as soon as 4 consecutive bases have lower mean quality'
    @params['minAvgQuality'] = 10
    @params['minReadLength'] = 19
    @params['name'] = 'CountSpacer_Result'
    @params['dictPath'] = ''
    @params['dictPath'] = {'select'=>''}
    Dir["/srv/GT/databases/GEML/sgRNA_Libs/*"].sort.select{|lib| File.directory?(lib)}.each do |dir|
      @params['dictPath'][File.basename(dir)] = File.basename(dir)
    end
    @params['leftPattern '] = ''
    @params['rightPattern'] = ''
    @params['maxMismatch'] = 1
    @params['maxMismatch', 'description'] = 'number of allowed mismatches for pattern search'

    @params['specialOptions'] = ''
    @params['specialOptions', 'description'] = 'special unsupported options that the R wrapper may support, format: <key>=<value>'
    @params['mail'] = ""
    @modules = ["QC/Flexbar", "QC/Trimmomatic", "Dev/Python2", "Dev/R", "Aligner/Bowtie"]
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
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@params['name'],
     'Report [File]'=>report_file,
     'Html [Link]'=>report_link,
     'Species'=>@dataset['Species'],
     'Read Count'=>@dataset['Read Count']
    }
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
