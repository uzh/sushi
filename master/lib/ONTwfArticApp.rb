#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class ONTwfArticApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'ONTwfArtic'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Variants'
    @description =<<-EOS
wf-artic: a Nextflow workflow for running the ARTIC SARS-CoV-2 workflow on multiplexed MinION, GridION, and PromethION runs<br/>
<a href='https://github.com/epi2me-labs/wf-artic'>ARTIC SARS-CoV-2 Workflow</a>
EOS

    @required_columns = ['Name', 'Read1', 'SampleSheet']
    @required_params = []
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '100'
    @params['normalise'] = '200'
    @params['normalise', 'description'] = 'depth ceiling for depth of coverage normalization'
    @params['schemeVersion'] = ['ARTIC/V3','ARTIC/V2','ARTIC/V1','ARTIC/V4','ARTIC/V4.1', 'Midnight-IDT/V1' ,'Midnight-ONT/V1','Midnight-ONT/V2','Midnight-ONT/V3','NEB-VarSkip/v1a-long','NEB-VarSkip/v1a','NEB-VarSkip/v2','NEB-VarSkip/v2b'] 
    @params['schemeVersion', 'description'] = 'Primer scheme version'
    @params['cmdOptions'] = ""
    @params['mail'] = ""
    @modules = ["Dev/jdk", "Dev/R"]
  end
  def next_dataset
       nds = {'Name'=>@dataset['Name']}     
       nds['ResultDir [File]'] = File.join(@result_dir, 'output/')
       nds['Report [Link]'] = File.join(@result_dir, 'output/wf-artic-report.html')
       nds
  end
  def commands
    run_RApp("EzAppONTwfArtic")
  end
end

if __FILE__ == $0

end

