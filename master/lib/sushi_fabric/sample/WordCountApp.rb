#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'

class WordCountApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Word_Count'
    @analysis_category = 'Stats'
    @required_columns = ['Name', 'Read1']
    @required_params = []
  end
  def next_dataset
    {'Name'=>@dataset['Name'],'Stats [File]'=>File.join(@result_dir, @dataset['Name'].to_s + '.stats')}
  end
  def preprocess
    @factors = get_columns_with_tag 'Factor'
    @factor_cols = @factors.first.keys
  end
  def commands
    commands = ''
    commands << "gunzip -c $GSTORE_DIR/#{@dataset['Read1']} |wc > #{@dataset['Name']}.stats\n"
    commands << "echo 'Factor columns: [#{@factor_cols.join(',')}]'\n"
    commands << "echo 'Factors: [#{@factors.join(',')}]'\n"
    commands
  end
end
if __FILE__ == $0
  usecase = WordCountApp.new

  usecase.project = "p1001"
  usecase.user = 'sushi_lover'
  usecase.parameterset_tsv_file = 'sample_parameterset.tsv'
  usecase.dataset_tsv_file = 'sample_dataset.tsv'
  #usecase.dataset_sushi_id = 26

  # run (submit to workflow_manager)
  usecase.run
  #usecase.test_run
end
