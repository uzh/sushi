#!/usr/bin/env ruby
# encoding: utf-8
Version = '20150617-032557'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class HTSeqApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'HTSeq'
    @analysis_category = 'Count'
    @required_columns = ['Name','BAM','BAI', 'refBuild']
    @required_params = ['build','paired', 'strandMode']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '16'
    @params['scratch'] = '100'
    @params['refBuild'] = ref_selector
    @params['paired'] = false
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['refFeatureFile'] = 'genes.gtf'
    @params['featureLevel'] = 'gene'
    @params['cmdOptions'] = ''
    @params['specialOptions'] = ''
    @params['mail'] = ""
  end
  def set_default_parameters
    @params['build'] = @dataset[0]['build']
    if dataset_has_column?('refFeatureFile')
      @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    end
    if dataset_has_column?('paired')
      @params['paired'] = @dataset[0]['paired']
    end                               
    if dataset_has_column?('strandMode')
      @params['strandMode'] = @dataset[0]['strandMode']
    end                               
  end

  def next_dataset
    {'Name'=>@dataset['Name'], 
     'Count [File]'=>File.join(@result_dir, "#{@dataset['Name']}.txt"), 
     'Species'=>@dataset['Species'],
     'build'=>@params['build'],
     'featureLevel'=>@params['featureLevel'],
     'refFeatureFile'=>@params['refFeatureFile'],
     'strandMode'=>@params['strandMode'],
     'paired'=>@params['paired'],
     'Read Count'=>@dataset['Read Count']
    }.merge(extract_column("Factor")).merge(extract_column("B-Fabric"))
  end
  def commands
    run_RApp("countHTSeqApp")
  end
end

if __FILE__ == $0
  usecase = HTSeqApp.new

  usecase.project = "p1001"
  usecase.user = 'masamasa'

  # set user parameter
  # for GUI sushi
  #usecase.params['process_mode'].value = 'SAMPLE'
  usecase.params['build'] = 'mm10'
  usecase.params['paired'] = true
  usecase.params['strandMode'] = 'both'
  usecase.params['cores'] = 8
  usecase.params['node'] = 'fgcz-c-048'

  # also possible to load a parameterset csv file
  # mainly for CUI sushi
  #usecase.parameterset_tsv_file = 'tophat_parameterset.tsv'
  #usecase.parameterset_tsv_file = 'test.tsv'

  # set input dataset
  # mainly for CUI sushi
  #usecase.dataset_tsv_file = 'tophat_dataset.tsv'

  # also possible to load a input dataset from Sushi DB
  usecase.dataset_sushi_id = 3

  # run (submit to workflow_manager)
  usecase.run
  #usecase.test_run

end

