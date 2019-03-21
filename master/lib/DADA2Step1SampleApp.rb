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
@params['paired'] = true
@params['paired', 'description'] = 'whether the reads are paired end; if false then only Read1 is considered even if Read2 is available.'
@params['mail'] = ""
@inherit_tags = ['B-Fabric', 'Characteristic', 'Mock','Group']
@modules = ['Dev/R']
end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
      if @params['Group']
      @required_columns << 'Group'
    end
        if @params['mockSample']
      @required_columns << 'mockSample'
    end
  end
 def set_default_parameters
    @params['paired'] = dataset_has_column?('Read2')
       @params['Group'] = dataset_has_column?('Group')
           @params['mockSample'] = dataset_has_column?('mockSample')
  end
  
def next_dataset
     nds = {'Name'=>@dataset['Name']}
     nds['RObjectWithSeqTab [File]'] = File.join(@result_dir, "#{@dataset['Name']}.seqTab.Rdata")
     nds['Technology [Factor]'] = @params['technology']
    pds = @dataset.clone
    pds.delete("Read1")
    pds.delete("Read2")
    pds.delete("Technology")
    nds.merge!(pds)
    nds
end
def commands
run_RApp("EzAppDADA2Step1Sample")
end
end

if __FILE__ == $0
end
