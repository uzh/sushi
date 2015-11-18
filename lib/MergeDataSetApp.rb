#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class MergeDataSetApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Merge_DataSet'
    @description = ""
    @description =<<-EOS
Merging two DataSets<br /><br />
For the moment, assuming that it works only on STARApp result datasets.<br />
<br />
Note
<ol>
<li>The number of samples between the two dataset should be same (not checked in the code).</li>
<li>The sample name should be same, namely the source dataset of STARApp should be same.</li>
<li>Actually, this is made for the input dataset of ReadClassifyApp.</li>
</ol>
    EOS
    @analysis_category = 'Demo'
    @params['DataSet'] = []
    @required_columns = ['Name', 'BAM', 'refBuild', 'Species']
    @required_params = ['DataSet']
  end
  def next_dataset
    {
      'Name'=>@dataset['Name'],
      'BAM1 [File]'=>@dataset['BAM'],
      'BAM2 [File]'=>@dataset['BAM2'],
      'refBuild1'=>@dataset['refBuild'],
      'refBuild2'=>@dataset['refBuild2'],
      'Species'=>@dataset['Species']
    }
  end
  def set_default_parameters
    if data_set = DataSet.find_by_id(@dataset_sushi_id)
      @params['DataSet'] = Hash[*data_set.project.data_sets.map{|d| [d.name, d.id]}.flatten]
    end
  end
  def preprocess
    dataset_sushi_id = @params['DataSet']
    dataset = DataSet.find_by_id(dataset_sushi_id)
    dataset_hash2 = {}
    if dataset = DataSet.find_by_id(dataset_sushi_id.to_i)
      dataset.samples.each do |sample|
        name = sample.to_hash['Name']
        dataset_hash2[name] = sample.to_hash
      end
    end
    # here
    # merge dataset_hash2 with @dataset_hash
    dataset_hash1 = @dataset_hash.clone
    dataset_hash1.each_with_index do |sample, i|
      name = sample['Name']
      @dataset_hash[i]['BAM2'] = dataset_hash2[name]['BAM [File]']
      @dataset_hash[i]['refBuild2'] = dataset_hash2[name]['refBuild']
    end
  end
  def commands
    "echo '#{GlobalVariables::SUSHI}'\n"
  end
end
