#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-095604'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class PostSamsa2AnalysisApp < SushiFabric::SushiApp
def initialize
  super
  @name = 'PostSamsa2Analysis'
  @analysis_category = 'Metagenomics'
  @description =<<-EOS
  Step_6 of Samsa2, data analysis and report. 
EOS
  @params['process_mode'] = 'DATASET'
  @required_columns = ['Name', 'annotationORGFileRefSeq', 'annotationFUNCFileRefSeq']
  @required_params = ['grouping', 'sampleGroup', 'refGroup'  ]
  @params['cores'] = '1'
  @params['cores', "context"] = "slurm"
  @params['ram'] = '7'
  @params['ram', "context"] = "slurm"
  @params['scratch'] = '10'
  @params['scratch', "context"] = "slurm"
  @params['grouping'] = '' 
  @params['grouping', "context"] = "PostSamsa2Analysis"
  @params['sampleGroup'] = '' 
  @params['sampleGroup', 'description'] = 'sampleGroup should be different from refGroup'
  @params['sampleGroup', "context"] = "PostSamsa2Analysis"
  @params['refGroup'] = '' 
  @params['refGroup', 'description'] = 'refGroup should be different from sampleGroup'
  @params['refGroup', "context"] = "PostSamsa2Analysis"
  @params['mail'] = ""
  @modules = ["Dev/R"]
    @inherit_columns = ["Order Id"]
end

def next_dataset
@params['name'] = "PostSamsa2Analysis"
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, '00index.html')
{'Name'=>@params['name'],
  'Report [File]'=>report_file,
  'Static Report [Link]'=>report_link,
}.merge(extract_columns(colnames: @inherit_columns))
end
def commands
run_RApp("EzAppPostSamsa2Analysis")
end
end

if __FILE__ == $0
end

