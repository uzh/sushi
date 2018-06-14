#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class MergeSCDataSetsApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'MergeSCDS'
    @description =<<-EOS
Merging more than two DataSets generaged from same library for CellRangerApp<br /><br />
Assuming that all other columns than file path are same between datasets.<br />
<br />
    EOS
    @analysis_category = 'SingleCell'
    @params['BaseDataSet'] = ''
    @params['TargetDataSet'] = []
    @required_columns = ['Name', 'Read1', 'RawDataDir', 'Species']
    @required_params = ['TargetDataSet']
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def next_dataset
    {
      'Name'=>@dataset['Name'],
      'Read1 [File]'=>@dataset['Read1'],
      'RawDataDir1 [File]'=>@dataset['RawDataDir'],
      'Read2 [File]'=>@dataset['Read1'],
      'RawDataDir2 [File]'=>@dataset['RawDataDir'],
      'Species'=>@dataset['Species'],
    }.merge(extract_columns(@inherit_tags))
  end
  def set_default_parameters
    if data_set = DataSet.find_by_id(@dataset_sushi_id)
      @params['BaseDataSet'] = data_set.name
      @params['TargetDataSet'] = Hash[*data_set.project.data_sets.map{|d| [d.name, d.id]}.flatten].to_a.reverse
    end
  end
  def preprocess
    dataset_sushi_id = @params['TargetDataSet']
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
    final_read_number = 1
    dataset_hash1.first.keys.each do |colname|
      if colname =~ /Read(\d+)\s+\[File\]/
        final_read_number = $1.to_i
      end
    end
    final_read_number += 1
    dataset_hash1.each_with_index do |sample, i|
      name = sample['Name']
      @dataset_hash[i]["Read#{final_read_number} [File]"] = dataset_hash2[name]['Read1 [File]']
      @dataset_hash[i]["RawDataDir#{final_read_number} [File]"] = dataset_hash2[name]['RawDataDir [File]']
    end
    @dataset_hash.sort_by!{|row| row['Name']}
  end
  def commands
    coms = ""
    coms << "echo '#{GlobalVariables::SUSHI}'\n"
    coms << "echo '#{GlobalVariables::SUSHI}' > #{@dataset['Name']}_dummy.txt\n"
    coms
  end
end
