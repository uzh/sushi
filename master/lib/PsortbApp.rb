#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-095008'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class PsortbApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Psortb'
    @analysis_category = 'Annotate'
    @description =<<-EOS
Psortb: subcellular protein localization prediction tool for Bacteria and Archea
<a href='https://www.psort.org/psortb/'>https://www.psort.org/psortb/</a>
EOS
    @required_columns = ['Name','Proteins']
    @required_params = ['cores', 'ram', 'scratch']
    # optional params
    @params['cores'] = '4'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '30'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '50'
    @params['scratch', "context"] = "slurm"
    @params['org'] = ['--negative', '--positive', '--archea']
    @params['org', 'description'] = 'type of organism: gram negative/ gram positive bacteria or archea'
    @params['org', "context"] = "Psortb"
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'specify other commandline options; do not specify any option that is already covered by the dedicated input fields'
    @params['cmdOptions', "context"] = "Psortb"
    @params['mail'] = ""
    @modules = ["Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'PsortbOut [File]'=>File.join(@result_dir, "#{@dataset['Name']}.psortb.txt"),
     'Proteins [File]'=>File.join(@result_dir, "#{@dataset['Name']}.proteins"),
     'Species'=>@dataset['Species'],
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppPsortb")
  end
end

if __FILE__ == $0

end

