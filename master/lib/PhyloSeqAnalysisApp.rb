#!/usr/bin/env ruby
# encoding: utf-8
Version = '20191108-134725'

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
@required_columns = ['Name', 'RObjectPhyloseq','RObjectQCChimera']
@required_params = ['taxonomicRank','rawCount','sampleFraction','numTopRanks']
@params['cores'] = '1'
@params['ram'] = '7'
@params['scratch'] = '10'
@params['taxonomicRank'] = ['Phylum','Class','Order','Family','Genus']
@params['taxonomicRank', 'description'] = 'Which rank should be displayed in the rank-dependent plots in the report?'
@params['rawCount'] = '5'
@params['rawCount', 'description'] = 'OTUs with fewer than these counts in less than the fraction of samples below will be removed.'
@params['sampleFraction'] = '0.3'
@params['sampleFraction', 'description'] = 'Minimum fraction of samples for which the raw count threshold above applies.'
@params['numTopRanks'] = '15'
@params['numTopRanks', 'description'] = 'Number of top-ranked OTUs to be considered in the plots.'
@params['grouping'] = ''
@params['grouping','description'] = 'This IS a selector. This decides from which of the lists below sampleGroup and refGroup must be chosen.'
@params['groupingDetails'] = ''
@params['groupingDetails','description'] = "This is NOT a selector. It simply shows the list of available values for the availalbe attributes  listed in grouping."
@params['sampleGroup'] = ''
@params['sampleGroup','description'] = 'This MUST be assigned from the values associated to the grouping variable'
@params['refGroup'] = ''
@params['refGroup','description'] = 'This MUST be assigned from the values associated to the grouping variable'
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
      ff = Hash[designMatrixTable.headers.collect{ |item| [item, designMatrixTable[item].uniq]}]
      @params['grouping']= ff.keys
      @params['groupingDetails']  = ff.to_a
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
