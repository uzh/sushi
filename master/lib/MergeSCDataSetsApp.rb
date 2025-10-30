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
    @params['FirstDataSet'] = ''
    @params['FirstDataSet', "context"] = "MergeSCDataSets"
    @params['SecondDataSet'] = ''
    @params['SecondDataSet', "context"] = "MergeSCDataSets"
    @required_columns = ['Name', 'Species', 'RawDataDir']
    @required_params = ['SecondDataSet']
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
    @child = true # child flag: true means that the next dataset is a child dataset
  end
  def next_dataset
    next_dataset_base = {
      'Name'=>@dataset['Name'],
    }
    @dataset.keys.select{|colname| colname =~ /RawDataDir\d*/ or colname =~ /Read\d+/}.each do |colname|
      new_colname = if colname == "RawDataDir"
                      "RawDataDir [File]"
                    else
                      "#{colname}"
                    end
      next_dataset_base[new_colname] = @dataset[colname]
    end
    next_dataset_base['Species'] = @dataset['Species']
    
    # Include Read Count column if it exists
    if @dataset['Read Count']
      next_dataset_base['Read Count'] = @dataset['Read Count']
    end
    
    next_dataset_base.merge(extract_columns(@inherit_tags))
  end

  def set_default_parameters
    if data_set = DataSet.find_by_id(@dataset_sushi_id)
      @params['FirstDataSet'] = data_set.name
      @params['SecondDataSet'] = Hash[*data_set.project.data_sets.map{|d| [d.name, d.id]}.flatten].to_a.reverse
    end
  end
  def preprocess
    dataset_sushi_id = @params['SecondDataSet']
    dataset = DataSet.find_by_id(dataset_sushi_id)
    dataset_hash2 = {}
    if dataset = DataSet.find_by_id(dataset_sushi_id.to_i)
      dataset.samples.each do |sample|
        name = sample.to_hash['Name']
        dataset_hash2[name] = sample.to_hash
      end
    end
    # here
    # merge dataset_hash2 with @dataset_hash
    dataset_hash1 = @dataset_hash.clone
    dataset_hash1.each_with_index do |sample, i|
      name = sample['Name']
      if dataset_hash2[name] and dataset_hash2[name]['RawDataDir [File]']
        @dataset_hash[i]["RawDataDir [File]"] += ",#{dataset_hash2[name]['RawDataDir [File]']}"
      end
      
      # Merge Read Count column by summing the values
      if dataset_hash2[name] and dataset_hash2[name]['Read Count']
        current_count = @dataset_hash[i]['Read Count'].to_i
        additional_count = dataset_hash2[name]['Read Count'].to_i
        @dataset_hash[i]['Read Count'] = (current_count + additional_count).to_s
      end
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
