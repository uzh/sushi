#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class EzPyzTestApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'EzPyzTestApp'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Testing'
    @description =<<-EOS
MultiSample Quality control after SpaceRanger<br/>
    EOS
    @required_columns = ['Name','IsTest']
    @required_params = ['name']
    @params['cores'] = '1'
    @params['ram'] = '30'
    @params['scratch'] = '100'
    @params['name'] = 'EzPyzTestApp'
    @params['sizeFactors'] = '1,3,5,10,20,30'
    @params['mail'] = ""
    @modules = ["Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def next_dataset
    report_dir = File.join(@result_dir, @params['name'])
    {'Name'=>@params['name'],
     'Report [File]'=>report_dir,
     'Report [Link]'=>File.join(report_dir, '00index.html'),
     'Species'=>(dataset = @dataset.first and dataset['Species'])
    }
  end
  def commands
    run_PyApp("EzPyzTestApp")
  end
end
