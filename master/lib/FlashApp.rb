#!/usr/bin/env ruby
# encoding: utf-8
Version = '20151104'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class FlashApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Flash'
    @analysis_category = 'Prep'
    @description =<<-EOS
Fast Length Adjustment of SHort reads
<a href='http://ccb.jhu.edu/software/FLASH/'>http://ccb.jhu.edu/software/FLASH/</a>
EOS
    
    @required_columns = ['Name','Read1','Read2']
    @required_params = ['paired']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '100'
    
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'specify the commandline options for flash; do not specify any option that is already covered by the dedicated input fields'
    @params['trimAdapter'] = true
    @params['trimAdapter', 'description'] = 'if adapters should be trimmed'
    @params['trimLeft'] = 0
    @params['trimLeft', 'description'] = 'fixed trimming at the "left" i.e. 5-prime end of the read'
    @params['trimRight'] = 0
    @params['trimRight', 'description'] = 'fixed trimming at the "right" i.e. 3-prime end of the read'
    @params['minTailQuality'] = 0
    @params['minTailQuality', 'description'] = 'if above zero, then reads are trimmed as soon as 4 consecutive bases have lower mean quality'
    @params['minAvgQuality'] = 10
    @params['minReadLength'] = 50
    @params['specialOptions'] = ''
    @params['specialOptions', 'description'] = 'special unsupported options that the R wrapper may support, format: <key>=<value>'
    @params['mail'] = ""
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
     'Read1 [File]'=>File.join(@result_dir, "#{@dataset['Name']}.extendedFrags.fastq.gz"), 
     'FlashLog [File]'=>File.join(@result_dir, "#{@dataset['Name']}_flash.log"),
     'TrimmomaticLog [File]'=>File.join(@result_dir, "#{@dataset['Name']}_preprocessing.log"),
     'Species'=>@dataset['Species'],
     'Read Count'=>@dataset['Read Count'],
     
     
    }.merge(extract_column("Factor")).merge(extract_column("B-Fabric"))
  end
  def commands
    run_RApp("EzAppFlash")
  end
end

if __FILE__ == $0
  run FlashApp
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

