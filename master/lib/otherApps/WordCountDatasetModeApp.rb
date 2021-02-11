#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class WordCountDatasetModeApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Word_Count_Dataset_Mode'
    @description = "test applicaiton #{GlobalVariables::SUSHI}"
    @analysis_category = 'Stats'
    @params['process_mode'] = 'DATASET'
    @params['cores'] = '1'
    @params['ram'] = '10'
    @params['scratch'] = '10'
    @params['custom_option'] = ['please select']
    @params['note'] = '' 
    @required_columns = ['Name', 'Read1']
    @required_params = []
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def next_dataset
    {
      'Name'=>@dataset[0]['Name'],
      'Stats [File]'=>File.join(@result_dir, @dataset[0]['Name'].to_s + '.stats')
    }.merge(extract_columns(@inherit_tags))
  end
  def set_default_parameters
    require 'csv'
    CSV.foreach("/srv/gstore/projects/p1001/ventricles/grouping.tsv", headers: true, col_sep: "\t") do |e|
      @params["custom_option"] << e["grouping"]
    end
  end
  def commands
    commands = ''
    commands << "echo '#{GlobalVariables::SUSHI}'\n"
    commands << "echo '#{SUSHI}'\n"
    commands
  end
end

if __FILE__ == $0
  run WordCountDatasetModeApp

  #usecase.project = "p1001"
  #usecase.user = 'sushi_lover'
  #usecase.parameterset_tsv_file = 'parameters.tsv'
  #usecase.dataset_tsv_file = 'dataset.tsv'
  #usecase.dataset_sushi_id = 1
  #usecase.params['grouping'] = 'hoge'

  # run (submit to workflow_manager)
  #usecase.run
  #usecase.test_run
end
