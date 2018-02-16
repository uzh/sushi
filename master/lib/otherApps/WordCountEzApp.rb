#!/usr/bin/env ruby
# encoding: utf-8
Version = '20180216-151217'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class WordCountEzApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'WordCountEz'
    @analysis_category = 'Test'
    @description =<<-EOS
Simplest SUSHIezRunApp
    EOS
    @required_columns = ['Name', 'Read1', 'Hoge']
    @required_params = ['refBuild']
    # optional params
    @params['cores'] = '1'
    @params['ram'] = '10'
    @params['scratch'] = '10'
    @params['refBuild'] = ref_selector
    @modules = ["Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'refBuild'=>@params['refBuild'],
     'Count [File]'=>File.join(@result_dir, "#{@dataset['Name']}_count.txt")
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppWordCount", "/srv/GT/analysis/masaomi/private_R_LIBS")
  end
end
