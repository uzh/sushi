#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-095604'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class MothurStep2DatasetApp < SushiFabric::SushiApp
def initialize
super
@name = 'MothurStep2Dataset'
@analysis_category = 'Metagenomics'
@description =<<-EOS
16S metagenomics data analysis and visualization with Phyloseq.
<a href='http://joey711.github.io/phyloseq/index.html'>http://joey711.github.io/phyloseq/index.html</a>
  EOS
@params['process_mode'] = 'DATASET'
@required_columns = ['Name', 'alignedFile', 'groupFile']
@required_params = ['cutOffTaxonomy', 'diffs','cutOffCluster']
@params['cutOffTaxonomy'] = '80'
@params['cutOffTaxonomy', 'description'] = 'Cut-off for taxonomy assignment in classify.seqs'
@params['diffs'] = '2'
@params['diffs', 'description'] = 'Differences allowed in the pre.cluster step. Should be 1 every 100 bases.'
@params['cutOffCluster'] = '0.03'
@params['cutOffCluster', 'description'] = 'Cut-off similarity to cluster OTUs'
@params['cores'] = '1'
@params['ram'] = '8'
@params['scratch'] = '10'
@params['representativeOTUs'] = ''
@params['representativeOTUs', 'description'] = 'Number of OTUs representing the sample.'
@params['mail'] = ""
@inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
@modules = ["Dev/R"]
end

  def preprocess
    if @params['Group']
      @required_columns << 'Group'
    end
  end
  def set_default_parameters
     @params['Group'] = dataset_has_column?('Group')
  end
  
def next_dataset
@params['name'] = "MothurStep2Dataset"
def next_dataset
     nds = {'Name'=>@dataset['Name']}
     nds['ChimeraPlot [File]'] = File.join(@result_dir, "#{@dataset['Name']}.chimPlot.txt")
     nds['PreClusteredAndChimeraSummary [File]'] = File.join(@result_dir, "#{@dataset['Name']}.preclChimSumm.txt")
      if @params['mockSample']
        nds['ErrorFile [File]'] = File.join(@result_dir, "#{@dataset['Name']}.errorCount.txt")
      end
     nds['stepConvergenceSummary [File]'] = File.join(@result_dir, "#{@dataset['Name']}.stepConv.txt")
     nds['OTUsToTaxonomyFile [File]'] = File.join(@result_dir, "#{@dataset['Name']}.OTUsToTax.txt")
     nds['OTUsCountTable [File]'] = File.join(@result_dir, "#{@dataset['Name']}.OTUsToCount.txt")
     nds

end
def commands
run_RApp("EzAppMothurStep2Dataset")
end
end

if __FILE__ == $0
end
