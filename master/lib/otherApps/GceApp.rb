#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class GceApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'GCE'
    @params['process_mode'] = 'SAMPLE'
    @analysis_category = 'QC'
    @description =<<-EOS
   GCE(genomic charactor estimator): a bayes model based method to estimate the genome size, genomic repeat content and the heterozygsis rate of the sequencing sample<br/>

    <a href='https://github.com/fanagislab/GCE'/>https://github.com/fanagislab/GCE</a> 
EOS
    @required_columns = ['Name','Read1']
    @required_params = ['cores']
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '100'
    @params['kSize'] = '17'
    @params['kSize', 'description'] = 'kmer size, recommand value 13 to 19, default=17'
    @params['cmdOptions'] = ""
    @params['mail'] = ""
    @modules = ["Dev/Perl/5.32.0", "Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'GCE Result [File]'=>File.join(@result_dir, "#{@dataset['Name']}"),
     'Report [Link]'=>File.join(@result_dir, "#{@dataset['Name']}", "#{@dataset['Name']}.html")
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppGce")
  end
end

if __FILE__ == $0
  
end
