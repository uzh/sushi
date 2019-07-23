#!/usr/bin/env ruby
# encoding: utf-8
Version = '20190723-101916'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class ConvERDSApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'ConvERDSApp'
    @description =<<-EOS
This converts the result DataSet of EAGLERCApp to an input DataSet of DNAHaplotypeCallerGVCFApp.
    EOS
    @analysis_category = 'Polyploid'
    @required_columns = ['Name', 'Species', 'dummy']
    @required_params = ['parent', 'type']
    # optional params
    @params['cores'] = '1'
    @params['ram'] = '10'
    @params['scratch'] = '10'
    @params['parent'] = ['1', '2']
    @params['type'] = ['ref', 'alt', 'unk']
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def preprocess
    dataset_hash = @dataset_hash.clone
    new_dataset_hash = []
    dataset_hash.each do |sample|
      new_sample_orig = {
        'Name'=>sample["Name"],
        'BAM'=>sample["Parent#{@params['parent']}RefBAM [File]"],
        'refBuild'=>sample["refBuild#{@params['parent']}"],
        'dummy'=>sample['dummy [File]'],
        'Species'=>sample['Species']
      }
      new_sample_other = {
        'Name'=>sample["Name"],
        'BAM'=>sample["Parent#{@params['parent']}AltBAM [File]"],
        'refBuild'=>sample["refBuild#{@params['parent']}"],
        'dummy'=>sample['dummy [File]'],
        'Species'=>sample['Species']
      }
      new_sample_common = {
        'Name'=>sample["Name"],
        'BAM'=>sample["Parent#{@params['parent']}UnkBAM [File]"],
        'refBuild'=>sample["refBuild#{@params['parent']}"],
        'dummy'=>sample['dummy [File]'],
        'Species'=>sample['Species']
      }

      if @params['type'] == 'ref'
        new_dataset_hash << new_sample_orig
      elsif @params['type'] == 'alt'
        new_dataset_hash << new_sample_other
      elsif @params['type'] == 'unk'
        new_dataset_hash << new_sample_common
      end
    end
    @dataset_hash = new_dataset_hash.sort_by{|row| row['Name']}
  end
  def next_dataset
    {'Name'=>@dataset['Name'], 
     'BAM'=>@dataset['BAM'],
     'refBuild'=>@dataset['refBuild'],
     'Species'=>@dataset['Species'],
     'Dummy [File]'=>File.join(@result_dir, "#{@dataset['Name']}_dummy.txt")
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    coms = ""
    coms << "echo '#{GlobalVariables::SUSHI}'\n"
    coms << "echo '#{GlobalVariables::SUSHI}' > #{@dataset['Name']}_dummy.txt\n"
    coms
  end
end


