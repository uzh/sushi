#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class EzPyzBin2CellApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'EzPyzBin2CellApp'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Spatial'
    @description =<<-EOS
    A Bin2Cell app for VisiumHD data<br/>
    EOS
    @required_columns = ['Name','BinnedOutputs2um','SourceImage','SpaceRanger']
    @required_params = ['name']
    @params['cores'] = '1'
    @params['ram'] = '30'
    @params['scratch'] = '100'
    @params['name'] = 'EzPyzBin2CellApp'
    @params['sizeFactors'] = '1,3,5,10,20,30'
    @params['mail'] = ""
    @modules = []
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def next_dataset
    report_dir = File.join(@result_dir, @params['name'])
    {'Name'=>@params['name'],
     'Figures [File]'=>File.join(report_dir, 'figures'),
     'Report [File]'=>File.join(report_dir, 'stardist')
    }
  end
  def commands
    run_PyApp("Bin2cell",conda_env: 'tmp_bin2cell')
  end
end
