#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-095008'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class SignalPApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'SignalP'
    @analysis_category = 'Annotate'
    @description =<<-EOS
SignalP: Signal peptide and cleavage sites in gram+, gram- and eukaryotic amino acid sequences
<a href='https://services.healthtech.dtu.dk/services/SignalP-5.0/'>https://services.healthtech.dtu.dk/services/SignalP-5.0/</a>
EOS
    @required_columns = ['Name','Proteins']
    @required_params = ['cores', 'ram', 'scratch']
    # optional params
    @params['cores'] = '4'
    @params['ram'] = '30'
    @params['scratch'] = '50'
    @params['org'] = ['gram-', 'gram+', 'arch', 'euk']
    @params['org', 'description'] = 'type of organism: gram negative/gram positive bacteria, Archaea, or Eukarya'
    @params['oformat'] = ['short', 'long']
    @params['oformat', 'description'] = 'Output format:long for generating the predictions with plots, short for the predictions without plots'
    @params['pformat'] = ['png', 'eps', 'none']
    @params['pformat', 'description'] = 'Plots output format. When long output selected, choose between png, eps or none to get just a tabular file'
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'specify other commandline options; do not specify any option that is already covered by the dedicated input fields'
    @params['mail'] = ""
    @modules = ["Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'Proteins [File]'=>File.join(@result_dir, "#{@dataset['Name']}.SignalP.faa"),
     'signalpOut [File]'=>File.join(@result_dir, "#{@dataset['Name']}.gff3"),
     'Species'=>@dataset['Species'],
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppSignalP")
  end
end

if __FILE__ == $0

end

