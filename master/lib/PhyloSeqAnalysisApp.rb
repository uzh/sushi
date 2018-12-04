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
@required_columns = ['Name', 'OTUsToTaxonomyFile', 'OTUsCountTable']
@required_params = ['representativeOTUs','group']
@params['cores'] = '1'
@params['ram'] = '8'
@params['scratch'] = '10'
@params['group'] = ['true','false']
@params['group', 'description'] = 'Is there a group factor in the dataset?'
@params['representativeOTUs'] = ''
@params['representativeOTUs', 'description'] = 'Number of OTUs representing the sample.'
@params['mail'] = ""
@inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
@modules = ["Dev/R"]
end
  
def next_dataset
@params['name'] = "PhyloSeq"
    report_file = File.join(@result_dir, "#{@params['name']}")
    report_link = File.join(report_file, '00index.html')
{'Name'=>@params['name'],
  'Report [File]'=>report_file,
  'Static Report [Link]'=>report_link,
}

end
def commands
run_RApp("EzAppPhyloSeqAnalysis")
end
end

if __FILE__ == $0
end
