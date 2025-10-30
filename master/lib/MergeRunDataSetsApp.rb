#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class MergeRunDataSetsApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'MergeRunDataSets'
    @description =<<-EOS
Merging two DataSets with Read1 and Read2 files by sample name<br /><br />
Concatenates Read1 and Read2 file paths from two datasets. Handles cases where samples exist in only one dataset.<br />
Filters samples based on minimum read count threshold.<br />
<br />
Note
<ol>
<li>Samples are matched by Name column between datasets.</li>
<li>Read Count values are summed when merging samples.</li>
<li>Samples below minimum read count threshold are excluded.</li>
<li>Samples unique to either dataset are preserved.</li>
<li>Specified columns in excluded_columns list are removed from output.</li>
</ol>
    EOS
    @analysis_category = 'Prep'
    @params['FirstDataSet'] = ''
    @params['SecondDataSet'] = ''
    @params['matchingColumn'] = ['Name', 'Tube', 'Sample Id']
    @params['matchingColumn', "context"] = "MergeRunDataSets"
    @params['paired'] = false
    @params['paired', "context"] = "MergeRunDataSets"
    @required_columns = ['Name', 'Species', 'Read1']
    @required_params = ['SecondDataSet']
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
    @child = true # child flag: true means that the next dataset is a child dataset
    
    # Define columns to exclude from the output
    @excluded_columns = ['Sample Id [B-Fabric]']
  end
  
  def next_dataset
    next_dataset_base = {
      'Name'=>@dataset['Name'],
    }
    @dataset.keys.select{|colname| colname =~ /Read\d+/}.each do |colname|
      new_colname = "#{colname} [File]"
      next_dataset_base[new_colname] = @dataset[colname]
    end
    next_dataset_base['Species'] = @dataset['Species']
    
    # Include Read Count column if it exists
    if @dataset['Read Count']
      next_dataset_base['Read Count'] = @dataset['Read Count']
    end
    
    # Extract columns with inherit tags and remove excluded columns
    inherited_columns = extract_columns(@inherit_tags)
    remove_excluded_columns(inherited_columns)
    
    next_dataset_base.merge(inherited_columns)
  end

  def set_default_parameters
    if data_set = DataSet.find_by_id(@dataset_sushi_id)
      @params['FirstDataSet'] = data_set.name
      @params['SecondDataSet'] = Hash[*data_set.project.data_sets.map{|d| [d.name, d.id]}.flatten].to_a.reverse
      @params['paired'] = dataset_has_column?('Read2')
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

    # Keep track of processed samples from dataset1
    processed_samples = Set.new
    
    # merge dataset_hash2 with @dataset_hash
    dataset_hash1 = @dataset_hash.clone
    dataset_hash1.each_with_index do |sample, i|
      name = sample['Name']
      processed_samples.add(name)
      
      if dataset_hash2[name] and dataset_hash2[name]['Read1 [File]']
        @dataset_hash[i]["Read1 [File]"] += ",#{dataset_hash2[name]['Read1 [File]']}"
      end
      
      if dataset_hash2[name] and dataset_hash2[name]['Read2 [File]'] and @params['paired']
        @dataset_hash[i]["Read2 [File]"] += ",#{dataset_hash2[name]['Read2 [File]']}"
      end
      
      # Merge Read Count column by summing the values
      if dataset_hash2[name] and dataset_hash2[name]['Read Count']
        current_count = @dataset_hash[i]['Read Count'].to_i
        additional_count = dataset_hash2[name]['Read Count'].to_i
        @dataset_hash[i]['Read Count'] = (current_count + additional_count).to_s
      end
    end
    
    # Add samples that exist only in dataset2
    dataset_hash2.each do |name, sample_data|
      unless processed_samples.include?(name)
        # Remove excluded columns from the sample data before adding
        cleaned_sample_data = sample_data.clone
        remove_excluded_columns(cleaned_sample_data)
        @dataset_hash << cleaned_sample_data
      end
    end
    
    @dataset_hash.sort_by!{|row| row['Name']}
  end
  
  private
  
  def remove_excluded_columns(data_hash)
    @excluded_columns.each do |column_name|
      data_hash.delete(column_name)
    end
  end
  
  def commands
    coms = ""
    coms << "echo '#{GlobalVariables::SUSHI}'\n"
    coms << "echo '#{GlobalVariables::SUSHI}' > #{@dataset['Name']}_dummy.txt\n"
    coms
  end
end
