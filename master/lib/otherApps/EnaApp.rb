#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class EnaApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'EnaApp'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'PublicData'
        @description =<<-EOS
Download public data from ENA<br/>
EOS
        @required_columns = ['Name','projectID']
        @required_params = ['name', 'projectID']
        @params['cores'] = '1'
        @params['ram'] = '15'
        @params['scratch'] = '300'
        @params['paired'] = false
        @params['name'] = 'ENA_Data'
        @params['projectID'] = ''
        @params['cmdOptions'] = ""
        @params['mail'] = ""
        @modules = ["Dev/R"]
  end
  def set_default_parameters
    @params['projectID'] = @dataset[0]['projectID']
  end
  def next_dataset
    {'Name'=>@dataset[0]['Name'],
     'projectID'=>@dataset[0]['projectID'],
     'ENA Result [File]'=>File.join(@result_dir, "#{@dataset[0]['Name']}"),
    }
  end
  def commands
    run_RApp("EzAppENA")
  end
end

if __FILE__ == $0
  usecase = FastqcApp.new

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

if __FILE__ == $0
  usecase = FastqcApp.new

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
