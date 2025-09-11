#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class EzPyzExampleApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'EzPyzExampleApp'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Testing'
    @description =<<-EOS
    A test app for ezPyz<br/>
    EOS
    @required_columns = ['Name','ResultDir','IsTest']
    @required_params = ['name']
    @params['cores'] = '1'
    @params['ram'] = '30'
    @params['scratch'] = '100'
    @params['name'] = 'EzPyzExampleApp'
    @params['sizeFactors'] = '1,3,5,10,20,30'
    @params['mail'] = ""
    @modules = []
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def next_dataset
    dir_name = "test_#{@params['name']}"
    report_dir = File.join(@result_dir, dir_name)
    {'Name'=> dir_name,
     'ReportData [File]'=>report_dir,
     'Report [Link]'=>File.join(report_dir, '00index.html'),
     'Species'=>(dataset = @dataset.first and dataset['Species'])
    }
  end
  def commands
    run_PyApp("Example",conda_env: 'tmp_ezpyz') #change to App specific env
  end
end
