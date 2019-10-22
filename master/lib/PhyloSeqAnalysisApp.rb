#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-095604'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class PhyloSeqAnalysisApp < SushiFabric::SushiApp
def initialize
super
@name = 'PhyloSeqAnalysis'
@analysis_category = 'Metagenomics'
@description =<<-EOS
16S metagenomics visualization with Phyloseq.
<a href='http://joey711.github.io/phyloseq/index.html'>http://joey711.github.io/phyloseq/index.html</a>
  EOS
@params['process_mode'] = 'DATASET'
@required_columns = ['Name', 'RObjectPhyloseq']
@required_params = ['representativeOTUs','group','taxonomicRank']
@params['cores'] = '1'
@params['ram'] = '8'
@params['scratch'] = '10'
@params['representativeOTUs'] = '150'
@params['representativeOTUs', 'description'] = 'Rough number of expected OTUs.'
@params['taxonomicRank'] = ['Phylum','Class','Order','Family','Genus']
@params['taxonomicRank', 'description'] = 'Which rank shuold be displayed in the rank-dependent plots in the report?'
@params['group'] = false
@params['group', 'description'] = 'are there group informations in the dataset?'
@params['rawCount'] = '5'
@params['rawCount', 'description'] = 'OTUs with fewer than these counts in less than the fraction of samples 
below will be removed.'
@params['sampleFraction'] = '0.3'
@params['sampleFraction', 'description'] = 'Minimum fraction of samples for which the raw count threshold above applies.'
@params['mail'] = ""
@params['name'] = "PhyloSeqReport"
@inherit_tags = ["B-Fabric", "Characteristic","group"]
@modules = ["Dev/R"]
end
def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, '00index.html')
    output_table = File.join(@result_dir, 'abundanceAndDiffExprTable.txt')
{'Name'=>@params['name'],
    'abundanceTable [File]'=>output_table,
    'Static Report [Link]'=>report_link,
    'Report [File]'=>report_file,
}
end
def commands
run_RApp("EzAppPhyloSeqAnalysis")
end
end

if __FILE__ == $0
end
