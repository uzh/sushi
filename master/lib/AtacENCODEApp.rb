#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class AtacENCODEApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'AtacENCODE'
    @params['process_mode'] = 'SAMPLE'
    @analysis_category = 'ATAC'
    @description =<<-EOS
    A ATAC-seq and DNase-seq processing pipeline from ENCODE. <br/>
    Fow now, it only supports human and mouse. <br/>

    <a href='https://github.com/ENCODE-DCC/atac-seq-pipeline'/>Github web-site</a>
    
EOS
    @required_columns = ['Name','Read1','Read2','Species']
    @required_params = ['name', 'paired']
    @params['cores'] = '16'
    @params['ram'] = '40'
    @params['scratch'] = '200'
    @params['paired'] = true
    @params['nReads'] = 25000000
    @params['name'] = 'AtacENCODE'
    @params['cmdOptions'] = ""
    @params['mail'] = ""
    @modules = ["Dev/jdk"]
  end
 def set_default_parameters
    @params['paired'] = dataset_has_column?('Read2')
  end
  def preprocess
    if @params['paired']
      @required_columns<<  'Read2'
    end
  end
  def next_dataset
     dataset = {
      'Name'=>@dataset['Name'],
      'Report [File]'=>File.join(@result_dir, "#{@dataset['Name']}_qc.html"),
      'Stats [File]'=>File.join(@result_dir, "#{@dataset['Name']}_qc.json"),
      'Html [Link]'=>File.join(@result_dir, "#{@dataset['Name']}_qc.html"),
    }
  end
  def commands
    run_RApp('EzAppAtacENCODE') #, conda_env: 'encode-atac-seq-pipeline')
  end
end

if __FILE__ == $0
  
end
