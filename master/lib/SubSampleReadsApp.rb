#!/usr/bin/env ruby
# encoding: utf-8
Version = '20220425'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class SubSampleReadsApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'SubSampleReads'
    @analysis_category = 'Prep'
    @description =<<-EOS
    R based subsampling of reads
EOS
    @required_columns = ['Name','Read1','Species','Pro']
    @required_params = ['paired', 'nReads']
    # optional params
    @params['cores'] = '1'
    @params['ram'] = '10'
    @params['scratch'] = '200'
    @params['paired'] = false
    @params['nReads'] = '500000'
    @params['mail'] = ""
    @modules = ["Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  end
  def set_default_parameters
    @params['paired'] = dataset_has_column?('Read2')
  end
  def next_dataset
     dataset = {
        'Name'=>@dataset['Name'],
        'Read1 [File]'=>File.join(@result_dir, "#{File.basename(@dataset['Read1'].to_s)}"),
        'Species'=>@dataset['Species'],
        'paired'=>@params['paired'],
        'Read Count'=>@dataset['Read Count']
    }.merge(extract_columns(@inherit_tags))
     if @params['paired'] 
       dataset['Read2 [File]'] = File.join(@result_dir, "#{File.basename(@dataset['Read2'].to_s)}")
     end
     dataset
  end
  def commands
    run_RApp("EzAppSubSampleReads")
  end
end

if __FILE__ == $0
  run SubSampleReadsApp

end
