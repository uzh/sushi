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
</ol>
    EOS
    @analysis_category = 'Prep'
    @params['FirstDataSet'] = ''
    @params['SecondDataSet'] = ''
    @params['matchingColumn'] = ['Name', 'Tube', 'Sample Id']
    @params['minReadCount'] = 10000
    @params['paired'] = false
    @required_columns = ['Name', 'Species', 'Read1']
    @required_params = ['SecondDataSet']
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
    @child = true # child flag: true means that the next dataset is a child dataset
  end
  
  def next_dataset
    next_dataset_base = {
      'Name'=>@dataset['Name'],
      'Species'=>@dataset['Species']
    }
    
    # Handle Read1 files (always present)
    next_dataset_base['Read1 [File]'] = @dataset['Read1 [File]']
    
    # Handle Read2 files if present and paired mode is enabled
    if @dataset['Read2 [File]'] && @params['paired']
      next_dataset_base['Read2 [File]'] = @dataset['Read2 [File]']
    end
    
    # Include Read Count if it exists
    if @dataset['Read Count']
      next_dataset_base['Read Count'] = @dataset['Read Count']
    end
    
    next_dataset_base.merge(extract_columns(@inherit_tags))
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
    
    # Get matching column (default to 'Name' if not found)
    match_col = @params['matchingColumn']
    if match_col.is_a?(Array)
      match_col = match_col.first
    end
    
    # Find samples that exist in both datasets
    dataset1_names = @dataset_hash.map { |sample| sample['Name'] }.compact
    dataset2_names = dataset_hash2.keys
    intersect_names = dataset1_names & dataset2_names
    unique_set1 = dataset1_names - dataset2_names
    unique_set2 = dataset2_names - dataset1_names
    
    puts "Files to merge: #{intersect_names.length}"
    puts "Files only in first dataset: #{unique_set1.length}"
    puts "Files only in second dataset: #{unique_set2.length}"
    
    # Process existing samples - merge with second dataset
    @dataset_hash.each_with_index do |sample, i|
      name = sample['Name']
      
      # Apply minimum read count filter
      read_count = sample['Read Count'].to_i
      if read_count < @params['minReadCount']
        @dataset_hash[i] = nil
        next
      end
      
      if intersect_names.include?(name) && dataset_hash2[name]
        # Merge Read1 files
        current_read1 = sample['Read1 [File]'] || ""
        additional_read1 = dataset_hash2[name]['Read1 [File]'] || ""
        
        if !additional_read1.empty?
          @dataset_hash[i]['Read1 [File]'] = current_read1.empty? ? additional_read1 : "#{current_read1},#{additional_read1}"
        end
        
        # Merge Read2 files if paired and both datasets have Read2
        if @params['paired'] && sample['Read2 [File]'] && dataset_hash2[name]['Read2 [File]']
          current_read2 = sample['Read2 [File]'] || ""
          additional_read2 = dataset_hash2[name]['Read2 [File]'] || ""
          
          if !additional_read2.empty?
            @dataset_hash[i]['Read2 [File]'] = current_read2.empty? ? additional_read2 : "#{current_read2},#{additional_read2}"
          end
        end
        
        # Sum Read Count values
        if dataset_hash2[name]['Read Count']
          current_count = sample['Read Count'].to_i
          additional_count = dataset_hash2[name]['Read Count'].to_i
          @dataset_hash[i]['Read Count'] = (current_count + additional_count).to_s
        end
      end
    end
    
    # Remove samples that didn't meet minimum read count
    @dataset_hash.compact!
    
    # Add samples that are unique to the second dataset
    unique_set2.each do |name|
      sample2 = dataset_hash2[name]
      
      # Apply minimum read count filter
      read_count = sample2['Read Count'].to_i
      next if read_count < @params['minReadCount']
      
      # Only add if Read2 requirement is met (if paired mode)
      if @params['paired'] && !sample2['Read2 [File]']
        puts "Skipping sample #{name} from second dataset: paired mode enabled but no Read2 file"
        next
      end
      
      @dataset_hash << sample2
    end
    
    # Handle the case where paired mode is requested but not all samples have Read2
    if @params['paired']
      samples_without_read2 = @dataset_hash.select { |sample| !sample['Read2 [File]'] || sample['Read2 [File]'].empty? }
      if !samples_without_read2.empty?
        puts "Warning: #{samples_without_read2.length} samples don't have Read2 files. Disabling paired mode."
        @params['paired'] = false
      end
    end
    
    @dataset_hash.sort_by! { |row| row['Name'] }
  end
  
  def commands
    coms = ""
    coms << "echo '#{GlobalVariables::SUSHI}'\n"
    coms << "echo '#{GlobalVariables::SUSHI}' > #{@dataset['Name']}_dummy.txt\n"
    coms
  end
end
