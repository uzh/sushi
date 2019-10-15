#!/usr/bin/env ruby
# encoding: utf-8
Version = '20180905-113400'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class DADA2Step1SampleApp < SushiFabric::SushiApp
def initialize
super
@name = 'DADA2Step1Sample'
@analysis_category = 'Metagenomics'
@description =<<-EOS
Data preprocssing with DADA2. Please make sure that the input files are from the same technology and adjust minLen and maxLen accordingly.
<a href='https://DADA2.org/wiki/MiSeq_SOP'>https://DADA2.org/wiki</a>
  EOS
  @params['process_mode'] = 'DATASET'
@required_columns = ['Name', 'Read1']
@required_params = ['maxLen', 'technology','paired']
@params['cores'] = '2'
@params['ram'] = '8'
@params['scratch'] = '10'
@params['maxLen'] = '300'
@params['maxLen', 'description'] = 'Sequences shorter than this long are removed.'
@params['technology'] = ['illumina','pacbio','ONT']
@params['technology', 'description'] = 'Sequencing technology used.'
@params['referenceFasta'] = ''
@params['referenceFasta', 'description'] = 'Full path to fasta file for the mock community (if available).'
@params['database'] = ['silva','RDP','greenGenes']
@params['database', 'description'] = '16S database to use for taxonomic assignment.'
@params['paired'] = true
@params['paired', 'description'] = 'whether the reads are paired end; if false then only Read1 is considered even if Read2 is available.'
@params['group'] = false
@params['group', 'description'] = 'is there at least one sample-group assignment column? If yes, ensure the 
    column name is in the format "NAME [Factor]"'
@params['mail'] = ""
@params['Name'] = "DADA2"
@inherit_tags = ['B-Fabric', 'Characteristic', 'group']
@modules = ['Dev/R']
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
      if @params['group']
     nds['sampleDescriptionFile [File]'] = File.join(@result_dir, "#{@params['Name']}.designMatrix.txt")
      end
     nds['OTUsToTaxonomyFile [File]'] = File.join(@result_dir, "#{@params['Name']}.OTUsToTax.txt")
     nds['OTUsCountTable [File]'] = File.join(@result_dir, "#{@params['Name']}.OTUsToCount.txt")
     nds['RObjectPhyloseq [File]'] = File.join(@result_dir, "#{@params['Name']}.phyloseq.Rdata")
     nds
end
def commands
run_RApp("EzAppDADA2Step1Sample")
end
end

if __FILE__ == $0
end
