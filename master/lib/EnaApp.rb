#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class EnaApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'EnaApp'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'PublicData'
        @description =<<-EOS
Download public data from ENA<br/>
EOS
        @required_columns = ['Name','projectID']
        @required_params = ['name', 'projectID']
        @params['cores'] = '1'
        @params['ram'] = '16'
        @params['scratch'] = '300'
        @params['paired'] = false
        @params['name'] = 'ENA_Data'
        @params['projectID'] = ''
        @params['cmdOptions'] = ""
        @params['mail'] = ""
        @modules = ["Dev/R"]
  end
  def set_default_parameters
    @params['projectID'] = @dataset[0]['projectID']
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'projectID'=>@dataset['projectID'],
     'ENA Result [File]'=>File.join(@result_dir, "#{@dataset['Name']}"),
    }
  end
  def commands
    run_RApp("EzAppENA")
  end
end

if __FILE__ == $0
  usecase = EnaApp.new

  usecase.project = "p1001"
  usecase.params['projectID'] = 'PRJEB12612'
  usecase.user = "lopitz"
end
