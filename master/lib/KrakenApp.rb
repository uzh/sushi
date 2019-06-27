#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class KrakenApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Kraken'
    @analysis_category = 'Metagenomics'
    @description =<<-EOS
Kraken taxonomic sequence classification system
<a href='https://ccb.jhu.edu/software/kraken2/index.shtml'>https://ccb.jhu.edu/software/kraken2/index.shtml</a>
EOS
    @required_columns = ['Name','Read1']
    @required_params = ['name', 'paired', 'krakenDBOpt']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '50'
    @params['trimAdapter'] = true
    @params['trimAdapter', 'description'] = 'if adapters should be trimmed'
    @params['trimLeft'] = 0
    @params['trimLeft', 'description'] = 'fixed trimming at the "left" i.e. 5-prime end of the read'
    @params['trimRight'] = 0
    @params['trimRight', 'description'] = 'fixed trimming at the "right" i.e. 3-prime end of the read'
    @params['minTailQuality'] = 20
    @params['minTailQuality', 'description'] = 'if above zero, then reads are trimmed as soon as 4 consecutive bases have lower mean quality'
    @params['minAvgQuality'] = 20
    @params['minReadLength'] = 50
    @params['krakenDBOpt'] = 'bacteria'
    @params['krakenDBOpt', 'description'] = 'kraken database options: viruses bacteria. Default is bac'
    @params['krakenConfidenceOpt'] = '0.0'
    @params['krakenConfidenceOpt', 'description'] = 'Confidence score threshold, between 0 and 1'
    @params['krakenPhredOpt'] = '0'
    @params['krakenPhredOpt', 'description'] = 'Phred score threshold'
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'specify other commandline options for kraken; do not specify any option that is already covered by the dedicated input fields'
    @params['mail'] = ""
    @modules = ["QC/Flexbar", "QC/Trimmomatic", "Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'KronaReport [Link]'=>File.join(@result_dir, "#{@dataset['Name']}.html"),
     'KrakenOut [File]'=>File.join(@result_dir, "#{@dataset['Name']}.txt"),
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppKraken")
  end
end

if __FILE__ == $0
end
