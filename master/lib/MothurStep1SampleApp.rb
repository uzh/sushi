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
@required_columns = ['Name', 'Read1']
@required_params = ['minLen', 'maxLen', 'technology','paired','mockSample']
@params['cores'] = '2'
@params['ram'] = '8'
@params['scratch'] = '10'
@params['minLen'] = '145'
@params['minLen', 'description'] = 'Sequences shorter than this long are removed.'
@params['maxLen'] = '330'
@params['maxLen', 'description'] = 'Sequences longer than this are removed.'
@params['technology'] = ['illumina','pacbio','ONT']
@params['technology', 'description'] = 'Sequencing technology used.'
@params['referenceFasta'] = ''
@params['referenceFasta', 'description'] = 'Full path to fasta file for the mock community (if available).'
@params['paired'] = false
@params['paired', 'description'] = 'whether the reads are paired end; if false then only Read1 is considered even if Read2 is available.'
@params['mail'] = ""
@inherit_tags = ['B-Fabric', 'Characteristic', 'Mock','Group']
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
      def preprocess
    if @params['Group']
      @required_columns << 'Group'
    end
  end
  def set_default_parameters
     @params['Group'] = dataset_has_column?('Group')
  end
def next_dataset
     nds = {'Name'=>@dataset['Name']}
     nds['RawDataSummary [File]'] = File.join(@result_dir, "#{@dataset['Name']}.rawSumm.txt")
     nds['DeduppedSummary [File]'] = File.join(@result_dir, "#{@dataset['Name']}.deduppedSumm.txt")
     nds['LenAndHomopSummary [File]'] = File.join(@result_dir, "#{@dataset['Name']}.lenHomopSumm.txt")
     nds['alignedFile [File]'] = File.join(@result_dir, "#{@dataset['Name']}.align.txt")
     nds['groupFile [File]'] = File.join(@result_dir, "#{@dataset['Name']}.group.txt")
     nds['Technology [Factor]'] = @params['technology']
    pds = @dataset.clone
    pds.delete("Read1")
    pds.delete("Read2")
    pds.delete("Technology")
    nds.merge!(pds)
    nds
end
def commands
run_RApp("EzAppMothurStep1Sample")
end
end

if __FILE__ == $0
end
