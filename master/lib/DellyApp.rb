#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class DellyApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'delly call'
    @analysis_category = 'Variants'
    @description =<<-EOS
delly uses paired-ends, split-reads and read-depth to sensitively and accurately delineate genomic rearrangements throughout the genome.<br/>
<a href='https://github.com/dellytools/delly'>Github website</a>
EOS

    @required_columns = ['Name','BAM','BAI', 'refBuild']
    @required_params = ['refBuild','ReadOpt']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '100'
    @params['refBuild'] = ref_selector
    @params['refBuild', 'description'] = 'the genome refBuild and annotation to use as reference.'
    @params['ReadOpt'] = 'sr'
    @params['ReadOpt', 'description'] = 'input read types: sr, pb, ont. Default is sr'
    @params['types'] = 'ALL'
    @params['types', 'description'] = 'SV types to call: DEL, INS, DUP, INV, BND, ALL. Default is ALL'
    @params['cmdOptions'] = ""
    @params['mail'] = ""
    @modules = ["Variants/SURVIVOR", "Dev/R"]
  end
  def next_dataset
    {
     'Name'=>@dataset['Name'],
     'refBuild'=>@params['refBuild'],
     'OutDir [File]'=>File.join(@result_dir, "#{@dataset['Name']}"),
     'OutReport [Link]'=>File.join(@result_dir, "#{@dataset['Name']}", "#{@dataset['Name']}.html")
    }
  end
  def set_default_parameters
      @params['refBuild'] = @dataset[0]['refBuild']
  end
  def commands
    run_RApp("EzAppDelly", conda_env: "tmp_delly")
  end
end

if __FILE__ == $0

end

