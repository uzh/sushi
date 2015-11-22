#!/usr/bin/env ruby
# encoding: utf-8
Version = '20151122-215452'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class ConvRCDSApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'ConvRCDSApp'
    @description =<<-EOS
This converts the result DataSet of ReadClassifyApp to an input DataSet of CountQCApp.
    EOS
    @analysis_category = 'Demo'
    @required_columns = ['Name','Species']
    @required_params = []
    # optional params
    @params['cores'] = '1'
    @params['ram'] = '10'
    @params['scratch'] = '10'
  end
  def preprocess
    dataset_hash = @dataset_hash.clone
    new_dataset_hash = []
    dataset_hash.each do |sample|
      new_sample_orig = {
        'Name'=>File.basename(sample['Parent1OrigBAM [File]']).gsub(/.bam/, ''),
        'BAM'=>sample['Parent1OrigBAM [File]'],
        'BAI'=>sample['Parent1OrigBAI [File]'],
        'refBuild'=>sample['refBuild1'],
        'Species'=>sample['Species']
      }
      new_sample_other = {
        'Name'=>File.basename(sample['Parent1OtherBAM [File]']).gsub(/.bam/, ''),
        'BAM'=>sample['Parent1OtherBAM [File]'],
        'BAI'=>sample['Parent1OtherBAI [File]'],
        'refBuild'=>sample['refBuild1'],
        'Species'=>sample['Species']
      }
      new_dataset_hash << new_sample_orig
      new_dataset_hash << new_sample_other
    end
    @dataset_hash = new_dataset_hash
  end
  def next_dataset
    {'Name'=>@dataset['Name'], 
     'BAM'=>@dataset['BAM'],
     'BAI'=>@dataset['BAI'],
     'refBuild'=>@dataset['refBuild'],
     'Species'=>@dataset['Species'],
     'dummy [File]'=>File.join(@result_dir, "#{@dataset['Name']}_dummy.txt")
    }.merge(extract_column("Factor")).merge(extract_column("B-Fabric"))
  end
  def commands
    coms = ""
    coms << "echo '#{GlobalVariables::SUSHI}'\n"
    coms << "echo '#{GlobalVariables::SUSHI}' > #{@dataset['Name']}_dummy.txt\n"
    coms
  end
end


