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
@required_params = ['group','taxonomicRank']
@params['cores'] = '1'
@params['ram'] = '8'
@params['scratch'] = '10'
@params['taxonomicRank'] = ['Phylum','Class','Order','Family','Genus']
@params['taxonomicRank', 'description'] = 'Which rank should be displayed in the rank-dependent plots in the report?'
@params['group'] = false
@params['group', 'description'] = 'are there group informations in the dataset?'
@params['rawCount'] = '5'
@params['rawCount', 'description'] = 'OTUs with fewer than these counts in less than the fraction of samples 
below will be removed.'
@params['sampleFraction'] = '0.3'
@params['sampleFraction', 'description'] = 'Minimum fraction of samples for which the raw count threshold above applies.'
@params['numTopRanks'] = '15'
@params['numTopRanks', 'description'] = 'Number of top-ranked OTUs to be considered in the plots.'
@params['mail'] = ""
@params['name'] = "PhyloSeqReport"
@inherit_tags = ["B-Fabric", "Characteristic","group"]
@modules = ["Dev/R"]
end
def set_default_parameters
      desMat=@dataset[0]['OTUsDesignMatrix [File]']
     if  !desMat.nil?
      require 'csv'
      fileName = File.join(SushiFabric::GSTORE_DIR, desMat)
      designMatrixTable = CSV.read(fileName,:headers => true, :col_sep => "\t")
      ff = Hash[designMatrixTable.headers.collect{ |item| [item, ['Values available for this variable']+designMatrixTable[item].uniq]}]
      ff.each do |key, value|
      @params[key] = value
      @params[key, 'description'] = 'This is NOT a selector, just a list of available values. If #{key} is used as grouping variable, sampleGroup and refGroup MUST be chose from this list.'
      end
      @params['grouping']= ff.keys
      @params['grouping','description']='This IS a selector. This decides from which of the lists above sampleGroup and refGroup must be chosen.'
      @params['sampleGroup'] = ''
      @params['sampleGroup','description'] = 'This MUST be chosen from the values in the dropdown list associated to the grouping variable'
      @params['refGroup'] = ''
      @params['refGroup','description'] = 'This MUST be a different choice from the values in the dropdown list associated to the SAME grouping variable'
      @required_params << ['sampleGroup','refGroup']
      end
end

def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, '00index.html')
{'Name'=>@params['name'],
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
