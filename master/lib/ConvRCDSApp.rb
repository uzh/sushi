#!/usr/bin/env ruby
# encoding: utf-8
Version = '20190723-101503'

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
    @analysis_category = 'Polyploid'
    @required_columns = ['Name','Species','dummy','Condition']
    @required_params = ['parent','type']
    # optional params
    @params['cores'] = '1'
    @params['ram'] = '10'
    @params['scratch'] = '10'
    @params['parent'] = ['1', '2']
    @params['type'] = ['all','w/o common','orig','other']
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def preprocess
    dataset_hash = @dataset_hash.clone
    new_dataset_hash = []
    dataset_hash.each do |sample|
      new_sample_orig = {
        'Name'=>File.basename(sample["Parent#{@params['parent']}OrigBAM [File]"]).gsub(/.bam/, ''),
        'BAM'=>sample["Parent#{@params['parent']}OrigBAM [File]"],
        'BAI'=>sample["Parent#{@params['parent']}OrigBAI [File]"],
        'refBuild'=>sample["refBuild#{@params['parent']}"],
        'Species'=>sample['Species'],
        'dummy'=>sample['dummy [File]'],
        'Condition'=>"#{sample['Condition [Factor]']}_orig"
      }
      new_sample_other = {
        'Name'=>File.basename(sample["Parent#{@params['parent']}OtherBAM [File]"]).gsub(/.bam/, ''),
        'BAM'=>sample["Parent#{@params['parent']}OtherBAM [File]"],
        'BAI'=>sample["Parent#{@params['parent']}OtherBAI [File]"],
        'refBuild'=>sample["refBuild#{@params['parent']}"],
        'Species'=>sample['Species'],
        'dummy'=>sample['dummy [File]'],
        'Condition'=>"#{sample['Condition [Factor]']}_other"
      }
      new_sample_common = {
        'Name'=>File.basename(sample["Parent#{@params['parent']}CommonBAM [File]"]).gsub(/.bam/, ''),
        'BAM'=>sample["Parent#{@params['parent']}CommonBAM [File]"],
        'BAI'=>sample["Parent#{@params['parent']}CommonBAI [File]"],
        'refBuild'=>sample["refBuild#{@params['parent']}"],
        'Species'=>sample['Species'],
        'dummy'=>sample['dummy [File]'],
        'Condition'=>"#{sample['Condition [Factor]']}_common"
      }

      unless @params['type'] == 'other'
        new_dataset_hash << new_sample_orig
      end
      unless @params['type'] == 'orig'
        new_dataset_hash << new_sample_other
      end
      if @params['type'] == 'all'
        new_dataset_hash << new_sample_common
      end
    end
    @dataset_hash = new_dataset_hash.sort_by{|row| row['Name']}
  end
  def next_dataset
    {'Name'=>@dataset['Name'], 
     'BAM'=>@dataset['BAM'],
     'BAI'=>@dataset['BAI'],
     'refBuild'=>@dataset['refBuild'],
     'Species'=>@dataset['Species'],
     'dummy [File]'=>File.join(@result_dir, "#{@dataset['Name']}_dummy.txt"),
     'Condition [Factor]'=>@dataset['Condition']
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    coms = ""
    coms << "echo '#{GlobalVariables::SUSHI}'\n"
    coms << "echo '#{GlobalVariables::SUSHI}' > #{@dataset['Name']}_dummy.txt\n"
    coms
  end
end


