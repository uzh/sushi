#!/usr/bin/env ruby
# encoding: utf-8
Version = '20180905-113400'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class MothurDataCleanBatchApp < SushiFabric::SushiApp
def initialize
super
@name = 'MothurDataCleanBatch'
@analysis_category = 'Metagenomics'
@description =<<-EOS
OTU-based metagenomics analysis with Mothur. Please make sure that the input files are from the same technology and adjust minLen and maxLen accordingly.
<a href='https://mothur.org/wiki/MiSeq_SOP'>https://mothur.org/wiki</a>
  EOS
@required_columns = ['Name', 'Read1']
@required_params = ['cutOffTaxonomy', 'diffs','cutOffCluster']
@params['cores'] = '1'
@params['ram'] = '8'
@params['scratch'] = '10'
@params['cutOffTaxonomy'] = '80'
@params['cutOffTaxonomy', 'description'] = 'Cut-off for taxonomy assignment in classify.seqs'
@params['diffs'] = '2'
@params['diffs', 'description'] = 'Differences allowed in the pre.cluster step. Should be 1 every 100 bases.'
@params['minLen'] = '145'
@params['minLen', 'description'] = 'Sequences shorter than this long are removed.'
@params['maxLen'] = '330'
@params['maxLen', 'description'] = 'Sequences longer than this are removed.'
@params['cutOffCluster'] = '0.03'
@params['cutOffCluster', 'description'] = 'Cut-off similarity to cluster OTUs'
@params['referenceFasta'] = '/srv/GT/analysis/grusso/courses/metagenomicsCourse/references/bioPoolReference.16S.fasta'
@params['referenceFasta', 'description'] = 'Full path to fasta file for the mock community (if available).'
@params['mail'] = ""
@inherit_tags = ["B-Fabric", "Characteristic", "Mock"]
@modules = ["Dev/R"]
end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  end
  def set_default_parameters
     @params['paired'] = dataset_has_column?('Read2')
  end
  def preprocess
    if @params['mockSample']
      @required_columns << 'Mock'
    end
  end
  def set_default_parameters
     @params['mockSample'] = dataset_has_column?('Mock')
  end
  
def next_dataset
     nds = {'Name'=>@dataset['Name']}
     nds['RawDataSummary [File]'] = File.join(@result_dir, "#{@dataset['Name']}.rawSumm.txt")
     nds['DeduppedSummary [File]'] = File.join(@result_dir, "#{@dataset['Name']}.deduppedSumm.txt")
     nds['LenAndHomopSummary [File]'] = File.join(@result_dir, "#{@dataset['Name']}.lenHomopSumm.txt")
     nds['MapFiltSummary [File]'] = File.join(@result_dir, "#{@dataset['Name']}.mapFilt.txt")
     nds['ChimeraPlot [File]'] = File.join(@result_dir, "#{@dataset['Name']}.chimPlot.txt")
     nds['PreClusteredAndChimeraSummary [File]'] = File.join(@result_dir, "#{@dataset['Name']}.preclChimSumm.txt")
      if @params['mockSample']
        nds['ErrorFile [File]'] = File.join(@result_dir, "#{@dataset['Name']}.errorCount.txt")
      end
     nds['stepConvergenceSummary [File]'] = File.join(@result_dir, "#{@dataset['Name']}.stepConv.txt")
     nds['OTUsToTaxonomyFile [File]'] = File.join(@result_dir, "#{@dataset['Name']}.OTUsToTax.txt")
     nds['OTUsCountTable [File]'] = File.join(@result_dir, "#{@dataset['Name']}.step.OTUsToCount.txt")
     nds = nds.merge(extract_columns(@inherit_tags))
end
def commands
run_RApp("EzAppMothurDataCleanBatch")
end
end

if __FILE__ == $0
end
