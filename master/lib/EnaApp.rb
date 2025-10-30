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
    @params['cores', \"context\"] = \"slurm\"
        @params['ram'] = '15'
    @params['ram', \"context\"] = \"slurm\"
        @params['scratch'] = '300'
    @params['scratch', \"context\"] = \"slurm\"
        @params['supportedMode'] = ['paired', 'single']
        @params['supportedMode', 'description'] = 'if mixed select either single read or paired end libraries from the study'
        @params['supportedMode', "context"] = "Ena"
        @params['tarOutput'] = false
        @params['tarOutput', 'description'] = 'for 10x data'
        @params['tarOutput', "context"] = "Ena"
        @params['excludedSamples'] = ''
        @params['excludedSamples', 'description'] = 'comma separated sample_accession SAMN-ids for samples to exclude'
        @params['excludedSamples', "context"] = "Ena"
        @params['includedSamples'] = ''
        @params['includedSamples', 'description'] = 'comma separated sample_accession SAMN-ids for samples to include'
        @params['includedSamples', "context"] = "Ena"
        @params['name'] = 'ENA_Data'
        @params['projectID'] = ''
        @params['projectID', "context"] = "Ena"
        @params['cmdOptions'] = ""
        @params['cmdOptions', "context"] = "Ena"
        @params['mail'] = ""
        @modules = ["Dev/R", "Dev/Python"]
  end
  def set_default_parameters
      @params['projectID'] = @dataset[0]['projectID']
      @params['datasetId'] = @dataset_sushi_id
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
