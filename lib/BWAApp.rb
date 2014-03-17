#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
require_relative 'optparse_ex'
include GlobalVariables

class BWAApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'BWA'
    @analysis_category = 'Map'
    @required_columns = ['Name','Read1','Species']
    @required_params = ['build','paired']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '16'
    @params['scratch'] = '100'
    @params['build'] = ref_selector
    @params['paired'] = false
    @params['algorithm'] = ['aln', 'mem', 'bwasw']
    @params['cmdOptions'] = ''
    @params['trimAdapter'] = false
    @params['trimLeft'] = 0
    @params['trimRight'] = 0
    @params['minTailQuality'] = 0
    @params['specialOptions'] = ''
  end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  end
  def set_default_parameters
    @params['paired'] = dataset_has_column?('Read2')
  end
  def next_dataset
    {'Name'=>@dataset['Name'], 
     'BAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam"), 
     'BAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam.bai"),
     'Species'=>@dataset['Species'],
     'build'=>@params['build'],
     'paired'=>@params['paired'],
     'Read Count'=>@dataset['Read Count']
    }.merge factor_dataset
  end
  def commands
    command = "/usr/local/ngseq/bin/R --vanilla --slave << EOT\n"
    command << "source('/usr/local/ngseq/opt/sushi_scripts/init.R')\n"
    command << "config = list()\n"
    config = @params
    config.keys.each do |key|
      command << "config[['#{key}']] = '#{config[key]}'\n" 
    end
    command << "config[['dataRoot']] = '#{@gstore_dir}'\n"
    command << "input = list()\n"
    input = @dataset
    input.keys.each do |key|
      command << "input[['#{key}']] = '#{input[key]}'\n" 
    end
    command << "output = list()\n"
    output = next_dataset
    output.keys.each do |key|
      command << "output[['#{key}']] = '#{output[key]}'\n" 
    end
    command << "mapBWA(input=input, output=output, config=config)\n"
    command << "EOT"
    command
  end
end

if __FILE__ == $0

end

