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
    @params['cores'] = '1'
    @params['ram'] = '4'
    @params['scratch'] = '100'

    @params['cmdOptions'] = '-no-g'
    @params['cmdOptions', 'description'] = 'specify the commandline options; do not specify any option that is already covered by the dedicated input fields'
    @params['dictPath'] = ''
    @params['keyStart'] = '24'
    @params['keyStart', 'description'] = 'start position of key'
    @params['keyEnd'] = '44'
    @params['keyEnd', 'description'] = 'end position of key'
    @params['specialOptions'] = ''
    @params['specialOptions', 'description'] = 'special unsupported options that the R wrapper may support, format: <key>=<value>'
    @params['mail'] = ""
    @modules = ["Dev/Python2", "Dev/R"]
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
    {'Name'=>@dataset['Name'],
     'Counts [File]'=>File.join(@result_dir, "#{@dataset['Name']}_counts.csv"),
     'Stats [File]'=>File.join(@result_dir, "#{@dataset['Name']}_statistics.txt")
     #'TrimmomaticLog [File]'=>File.join(@result_dir, "#{@dataset['Name']}_preprocessing.log"),
     #'Species'=>@dataset['Species'],
     #'Read Count'=>@dataset['Read Count'],
    }#.merge(extract_columns(@inherit_tags))
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
