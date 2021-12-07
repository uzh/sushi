#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

require 'csv'
class WordCountApp < SushiFabric::SushiApp
  def initialize
    super
    @employee = true
    @name = 'Word_Count'
    @description = "test applicaiton #{GlobalVariables::SUSHI}"
    @analysis_category = 'Stats'
    @params['cores'] = '1'
    @params['ram'] = '10'
    @params['scratch'] = '10'
    @params['grouping'] = '' ### TODO: this should be filled by a column selector that allows to select a column with the tag 'Factor'
    @params['sampleGroup'] = '' ## TODO: this should be a value from the selected column
    @params['refGroup'] = '' ## TODO: this should be a value from the selected column
    @params['count_option'] = ['', '-c', '-l', '-m', '-w']
    @params['test_option1'] = ['option1', 'option2', 'option3']
    @params['test_option1', 'multi_selection'] = true
    @params['test_option1', 'multi_selection_size'] = 10
    @params['test_option2'] = ''
    @params['test_option2', 'multi_selection'] = true
    @params['test_option2', 'all_selected'] = true
    @params['test_option3'] = ''
    @params['test_option3', 'file_upload'] = true
    @params['cellType'] = []
    @params['cellType', 'multi_selection'] = true
    @params['note'] = '' 
    @params['note', 'employee'] = true
    @required_columns = ['Name', 'Read1']
    @required_params = []
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def next_dataset
    {
      'Name'=>@dataset['Name'],
      'Stats [File]'=>File.join(@result_dir, @dataset['Name'].to_s + '.stats'),
      'Options' => @params['test_option1']
    }.merge(extract_columns(@inherit_tags))
  end
  def set_default_parameters
    if @dataset[0]['Options']
      @params['test_option2'] = @dataset[0]['Options'].split(',')
    end
    cell_types = {}
    CSV.foreach("/srv/GT/databases/all_cell_markers.txt", headers: true, col_sep: "\t") do |e|
      cell_types[e["cellType"]] = true
    end
    @params['cellType'] = cell_types.keys.sort
    @params['note'] = @dataset[0]['Read1 [File]']
  end
  def preprocess
    if @factors = get_columns_with_tag('Factor') and @factors.first
      @factor_cols = @factors.first.keys
    end
  end
  def commands
    commands = ''
    commands << "gunzip -c $GSTORE_DIR/#{@dataset['Read1']} |wc > #{@dataset['Name']}.stats\n"
    commands << "echo 'Factor columns: [#{@factor_cols.join(',')}]'\n"
    commands << "echo 'Factors: [#{@factors.join(',')}]'\n"
    commands << "echo 'test_option1: #{@params["test_option1"]}'\n"
    commands << "echo 'test_option2: #{@params["test_option2"]}'\n"
    commands << "echo '#{GlobalVariables::SUSHI}'\n"
    commands << "echo '#{SUSHI}'\n"
    commands
  end
end

if __FILE__ == $0
  run WordCountApp

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
