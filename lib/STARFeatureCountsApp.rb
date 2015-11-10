#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

require 'STARApp'
require 'FeatureCountsApp'

class STARFeatureCountsApp < SushiFabric::SushiApp
  def initialize
    super

    sushi_apps = [STARApp, FeatureCountsApp]
    @apps = sushi_apps.map{|klass| klass.new}
    @name = self.class.to_s
    @description = @apps.map{|app| app.description}.join("\n")
    @analysis_category = @apps.map{|app| app.analysis_category}.join

    @params = @apps.first.params
    @apps.last.params.each do |key, value|
      @params[key] = value
    end
    @apps.last.instance_variable_set(:@params, @params)
    @required_columns = @apps.first.required_columns
    @required_params = @apps.first.required_params
  end
  def preprocess
    @apps.first.dataset_sushi_id = @dataset_sushi_id
    @apps.first.project = @project
    @apps.first.instance_variable_set(:@name, @name)
    @apps.first.next_dataset_name = @next_dataset_name
    @apps.first.user = @user
    @apps.first.parameterset_tsv_file = @parameterset_tsv_file

    @apps.first.set_input_dataset
    @apps.first.set_dir_paths
    @apps.first.preprocess
    @apps.first.set_output_files
    @apps.first.set_user_parameters

    # merge next_dataset to input_dataset
    @apps.first.dataset_hash.each_with_index do |row, i|
      dataset = Hash[*row.map{|key,value| [key.gsub(/\[.+\]/,'').strip, value]}.flatten]
      @apps.first.instance_variable_set(:@dataset, dataset)
      @dataset_hash[i].merge!(@apps.first.next_dataset) 
    end

    #@apps.last.dataset_sushi_id = @dataset_sushi_id
    @apps.last.project = @project
    @apps.last.instance_variable_set(:@name, @name)
    @apps.last.next_dataset_name = @next_dataset_name
    @apps.last.user = @user
    @apps.last.parameterset_tsv_file = @parameterset_tsv_file

    #@apps.last.set_input_dataset
    @apps.last.instance_variable_set(:@dataset_hash, @dataset_hash)
    @apps.last.set_dir_paths
    @apps.last.preprocess
    @apps.last.set_output_files
    @apps.last.set_user_parameters
  end
  def next_dataset
#    @apps.first.instance_variable_set(:@dataset, @dataset)
#    @apps.first.instance_variable_set(:@result_dir, @result_dir)
#    @apps.first.next_dataset

    @apps.last.instance_variable_set(:@dataset, @dataset)
    @apps.last.instance_variable_set(:@result_dir, @result_dir)
    @apps.last.next_dataset
  end
  def jobfooter(output_files, next_dataset)
    footer = ""
    output_files.map{|header| next_dataset[header]}.each do |file|
      src_file = File.basename(file)
      dest_dir = File.dirname(File.join(@gstore_dir, file))
      footer << copy_commands(src_file, dest_dir).join("\n")+"\n"
    end
    footer
  end
  def commands
    @apps.first.instance_variable_set(:@dataset, @dataset)
    @apps.first.instance_variable_set(:@result_dir, @result_dir)
    @apps.last.instance_variable_set(:@dataset, @dataset)
    @apps.last.instance_variable_set(:@result_dir, @result_dir)
    
    [@apps.first.commands, jobfooter(@apps.first.instance_variable_get(:@output_files), @apps.first.next_dataset), @apps.last.commands].join("\n") 
  end
end

