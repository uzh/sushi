#!/usr/bin/env ruby
# encoding: utf-8
Version = '20181127-163024'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class MothurStep2DatasetReportApp < SushiFabric::SushiApp
def initialize
super
@name = 'MothurStep2DatasetReport'
@analysis_category = 'Metagenomics'
@description =<<-EOS
Report of the data analysis performed with Mothur.
EOS
@params['process_mode'] = 'DATASET'
@required_columns = ['Name','ChimeraPlot','PreClusteredAndChimeraSummary','stepConvergenceSummary','OTUsToTaxonomyFile','OTUsCountTable']
@required_params = ['name', 'group']
@params['cores'] = '1'
@params['ram'] = '8'
@params['scratch'] = '10'
@params['mail'] = ""
@params['Group'] = 'true'
@params['Group', 'description'] = 'Is there a group factor in the dataset?'
@params['name'] = "EzAppMothurStep2DatasetReport"
@inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
@modules = ["Dev/R/3.5.1"]
end
def next_dataset
report_file = File.join(@result_dir, '00index_files')
report_link = File.join(@result_dir, '00index.html')
{'Name'=>@params['name'],
  'Report [File]'=>report_file,
  'Static Report [Link,File]'=>report_link
}
end
def commands
run_RApp("EzAppMothurStep2DatasetReport")
end
end

if __FILE__ == $0
end
