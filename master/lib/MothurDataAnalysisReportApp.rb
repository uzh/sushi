#!/usr/bin/env ruby
# encoding: utf-8
Version = '20181127-160812'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class MothurDataAnalysisReportApp < SushiFabric::SushiApp
def initialize
super
@name = 'MothurDataAnalysisReport'
@analysis_category = 'Metagenomics'
@description =<<-EOS
Report of the data analysis performed with Mothur.
EOS
@params['process_mode'] = 'DATASET'
@required_columns = ['Name','RawDataSummary','DeduppedSummary','LenAndHomopSummary','MapFiltSummary','ChimeraPlot','PreClusteredAndChimeraSummary','stepConvergence']
@required_params = ['name']
@params['cores'] = '1'
@params['ram'] = '8'
@params['scratch'] = '10'
@params['mail'] = ""
@params['name'] = "MothurDataAnalysisReport"
@inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
@modules = ["Dev/R"]
end
def next_dataset
report_file = File.join(@result_dir, '00index_files')
report_link = File.join(@result_dir, '00index.html')
{'Name'=>@dataset['Name'],
  'Report [File]'=>report_file,
  'Static Report [Link,File]'=>report_link,
}

end
def commands
run_RApp("EzAppMothurDataAnalysisReport")
end
end

if __FILE__ == $0
end
