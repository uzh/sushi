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
@params['representativeOTUs'] = ''
@params['representativeOTUs', 'description'] = 'Number of OTUs representing the sample.'
@params['mail'] = ""
@params['name'] = "PhyloSeqReport"
@inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
@modules = ["Dev/R"]
end
      def preprocess
    if @params['group']
      @required_columns << 'sampleDescriptionFile'
    end
  end
  def set_default_parameters
     @params['group'] = dataset_has_column?('sampleDescriptionFile')
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
