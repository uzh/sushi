#!/usr/bin/env ruby
# encoding: utf-8
Version = '20181127-163024'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class MothurStep1SampleReportApp < SushiFabric::SushiApp
def initialize
super
@name = 'MothurStep1SampleReport'
@analysis_category = 'Metagenomics'
@description =<<-EOS
Report of the preprocessin analysis performed with Mothur.
EOS
@params['process_mode'] = 'DATASET'
@required_columns = ['Name','RawDataSummary','DeduppedSummary','LenAndHomopSummary']
@required_params = ['name']
@params['cores'] = '1'
@params['ram'] = '7'
@params['scratch'] = '10'
@params['mail'] = ""
@params['name'] = "MothurStep1SampleReportApp"
@inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
@modules = ["Dev/R/3.5.1"]
end
def next_dataset
@params['name'] = "MothurStep1"
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, '00index.html')
{'Name'=>@params['name'],
  'Report [File]'=>report_file,
  'Static Report [Link]'=>report_link,
}
end
def commands
run_RApp("EzAppMothurStep1SampleReport")
end
end

if __FILE__ == $0
end
