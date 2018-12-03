#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-095422'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class SpadesApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Spades'
    @analysis_category = 'Assemble'
    @description =<<-EOS
SPAdes genome assembler
<a href='http://cab.spbu.ru/software/spades/'>http://cab.spbu.ru/software/spades/</a>
EOS

    @required_columns = ['Name','Read1']
    @required_params = ['cores', 'ram', 'scratch']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '50'
    @params['scratch'] = '100'
    @params['paired'] = true
    @params['genomeType'] = ["isolate","metagenome"]
    @params['genomeType', 'description'] = 'are you assemblying an isolate or a metagenome?'
    @params['otherBasicOpt'] = ''
    @params['otherBasicOpt', 'description'] = 'SPAdes basic options (apart from --meta): --sc, --rna, --plasmid, Default is empty for genome assembly without MDA'
    @params['spadesPipeOpt'] = '--careful'
    @params['spadesPipeOpt', 'description'] = 'SPAdes pipeline options: --only-assembler, --careful'
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'specify other commandline options for SPAdes; do not specify any option that is already covered by the dedicated input fields'
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
    @params['mail'] = ""
    @modules = ["QC/Flexbar", "Assembly/SPAdes", "QC/Trimmomatic", "Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'Draft [File]'=>File.join(@result_dir, "#{@dataset['Name']}.fasta"),
     'SpadesOut [File]'=>File.join(@result_dir, "#{@dataset['Name']}"),
     'SpadesLog [File]'=>File.join(@result_dir, "#{@dataset['Name']}_spades.log"),
     'TrimmomaticLog [File]'=>File.join(@result_dir, "#{@dataset['Name']}_preprocessing.log"),
     'Species'=>@dataset['Species'],
     'Read Count'=>@dataset['Read Count'],
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppSpades")
  end
end

if __FILE__ == $0
end
