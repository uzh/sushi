#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class MageckTestApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'MageckTest'
    @analysis_category = 'QC'
    @description =<<-EOS
    Run test module in the tool Model-based Analysis of Genome-wide CRISPR-Cas9 Knockout (<a href='https://sourceforge.net/p/mageck/wiki/Home/'>MAGeCK</a>)
    EOS
    @params['process_mode'] = 'DATASET'
    @required_columns = ['Name','Count']
    @required_params = []
    # optional params
    @params['cores'] = ['1']
    @params['ram'] = ['4']
    @params['scratch'] = ['10']
    @params['name'] = 'MAGeCK_Test'
    @params['species'] = ['hsa', 'mmu']
    @params['specialOptions'] = ''
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'specify other commandline options for MAGeCK_Test; do not specify any option that is already covered by the dedicated input fields'
    @params['grouping'] = ''
    @params['sampleGroup'] = ''
    @params['sampleGroup', 'description'] = 'sampleGroup should be different from refGroup'
    @params['refGroup'] = ''
    @params['refGroup', 'description'] = 'refGroup should be different from sampleGroup'
    @params['mail'] = ""
    @modules = ["Dev/R"]
    @inherit_tags = ["B-Fabric"]
  end
  def preprocess
  end
  def next_dataset
    @comparison = "#{@params['sampleGroup']}--over--#{@params['refGroup']}"
    @params['comparison'] = @comparison
    @params['name'] = @comparison
    report_folder = File.join(@result_dir, "#{@params['comparison']}")
    {'Name'=>@params['name'],
     'Report [File]'=>report_folder,
    }.merge(extract_columns(@inherit_tags))
  end
  def set_default_parameters
  end
  def commands
    run_RApp("EzAppMageckTest")
  end
end

if __FILE__ == $0

end

