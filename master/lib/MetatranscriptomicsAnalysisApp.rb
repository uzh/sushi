#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-095604'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class MetatranscriptomicsAnalysisApp < SushiFabric::SushiApp
def initialize
super
@name = 'MetatranscriptomicsAnalysis'
@analysis_category = 'Metagenomics'
@description =<<-EOS
Analysis of metatranscriptomics data produced by Samsa2. 
  EOS
@params['process_mode'] = 'DATASET'
@required_columns = ['Name', 'annotationFileRefSeq']
@required_params = ['numberOfTopNCategories','grouping', 'sampleGroup', 'refGroup'  ]
@params['cores'] = '1'
@params['ram'] = '8'
@params['scratch'] = '10'
@params['grouping'] = '' 
@params['sampleGroup'] = '' 
@params['sampleGroup', 'description'] = 'sampleGroup should be different from refGroup'
@params['refGroup'] = '' 
@params['refGroup', 'description'] = 'refGroup should be different from sampleGroup'
@params['numberOfTopNCategories'] = '30'
@params['numberOfTopNCategories', 'description'] = 'Number of top most represented functions to report.'
@params['mail'] = ""
@inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
@modules = ["Dev/R"]
end

def next_dataset
@params['name'] = "MetatranscriptomicsAnalysis"
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, '00index.html')
{'Name'=>@params['name'],
  'Report [File]'=>report_file,
  'Static Report [Link]'=>report_link,
}
end
def commands
run_RApp("EzAppMetatranscriptomicsAnalysis")
end
end

if __FILE__ == $0
end
