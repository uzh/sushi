#!/usr/bin/env ruby
# encoding: utf-8
Version = '20180905-113400'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class MothurStep1SampleApp < SushiFabric::SushiApp
def initialize
super
@name = 'MothurStep1Sample'
@analysis_category = 'Metagenomics'
@description =<<-EOS
Data preprocssing with Mothur. Please make sure that the input files are from the same technology and adjust minLen and maxLen accordingly.
<a href='https://mothur.org/wiki/MiSeq_SOP'>https://mothur.org/wiki</a>
  EOS
   @params['process_mode'] = 'DATASET'
@required_columns = ['Name', 'Read1', 'Adapter1']
@required_params = ['technology','paired','cutOffTaxonomy','diffs','cutOffCluster']
@params['cores'] = '2'
@params['ram'] = '8'
@params['scratch'] = '10'
@params['technology'] = ['illumina','pacbio']
@params['technology', 'description'] = 'Sequencing technology used.'
@params['referenceFasta'] = ''
@params['referenceFasta', 'description'] = 'Full path to fasta file for the mock community (if available).'
@params['paired'] = false
@params['paired', 'description'] = 'whether the reads are paired end; if false then only Read1 is considered even if Read2 is available.'
@params['diffs'] = '2'
@params['diffs', 'description'] = 'Differences allowed in the pre.cluster step. Should be 1 every 100 bases.'
@params['cutOffTaxonomy'] = '80'
@params['cutOffTaxonomy', 'description'] = 'Cut-off for taxonomy assignment in classify.seqs'
@params['cutOffCluster'] = '0.03'
@params['cutOffCluster', 'description'] = 'Cut-off similarity to cluster OTUs'
@params['group'] = false
@params['group', 'description'] = 'is there at least one sample-group assignment column? If yes, ensure the 
    column name is in the format "NAME [Factor]"'
@params['mail'] = ""
@params['Name'] = "Mothur"
@inherit_tags = ['B-Fabric', 'Characteristic', 'Mock','Group']
@modules = ['Dev/R', 'Tools/Mothur']
end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  end
 def set_default_parameters
    @params['paired'] = dataset_has_column?('Read2')
  end
  
def next_dataset
     nds = {'Name'=>@params['Name']}
     nds['OTUsToTaxonomyFile [File]'] = File.join(@result_dir, "#{@params['Name']}.OTUsToTax.txt")
     nds['OTUsCountTable [File]'] = File.join(@result_dir, "#{@params['Name']}.OTUsToCount.txt")
     if @params['group']
     nds['OTUsDesignMatrix [File]'] = File.join(@result_dir, "#{@params['Name']}.designMatrix.txt")
     end
     nds['RObjectPhyloseq [File]'] = File.join(@result_dir, "#{@params['Name']}.phyloseq.Rdata")
     nds['RObjectQCChimera [File]'] = File.join(@result_dir, "#{@params['Name']}.QCChimera.Rdata")
     nds
end
def commands
run_RApp("EzAppMothurStep1Sample")
end
end

if __FILE__ == $0
end
