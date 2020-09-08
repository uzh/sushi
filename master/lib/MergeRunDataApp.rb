#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class MergeRunDataApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'MergeRunData'
    @description =<<-EOS
Merging fastq files from two illumina runs by name <br /><br />
<br />
    EOS
    @analysis_category = 'Prep'
    @params['process_mode'] = 'DATASET'
    @params['FirstDataSet'] = ''
    @params['SecondDataSet'] = ''
    @params['Name'] = 'MergedRunData'
    @required_columns = ['Name', 'Species', 'Read1']
    @required_params = ['SecondDataSet', 'paired']
    
    # optional params
    @params['cores'] = '1'
    @params['ram'] = '4'
    @params['scratch'] = '200'
    @params['paired'] = false
    @params['paired', 'description'] = 'either the reads are paired-ends or single-end'
    @params['mail'] = ''
    
    @modules = ["Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
    @child = false # child flag: true means that the next dataset is a child dataset
  end
  def next_dataset
    report_file = File.join(@result_dir, @params['Name'])
    {'Name'=>@params['Name'],
     'Result [File]'=>report_file,
    }
  end

  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  end
  
  def set_default_parameters
    if data_set = DataSet.find_by_id(@dataset_sushi_id)
      @params['FirstDataSet'] = data_set.name
      @params['SecondDataSet'] = Hash[*data_set.project.data_sets.map{|d| [d.name, d.id]}.flatten].to_a.reverse
      @params['paired'] = dataset_has_column?('Read2')
    end
  end
  
  def commands
     dataset2 = DataSet.find_by_id(params['SecondDataSet'])
     @params['DataSetName2'] = dataset2.name
     run_RApp("EzAppMergeRunData")
  end
end
