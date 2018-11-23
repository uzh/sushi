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
OTU-based metagenomics analysis with Mothur.
<a href='https://mothur.org/wiki/MiSeq_SOP'>https://mothur.org/wiki</a>
  EOS
@required_columns = ['Name', 'Read1']
@required_params = ['cutOff', 'diffs']
@params['cores'] = '1'
@params['ram'] = '8'
@params['scratch'] = '10'
@params['cutOff'] = '80'
@params['cutOff', 'description'] = 'Cut-off for taxonomy assignment in classify.seqs'
@params['diffs'] = '2'
@params['diffs', 'description'] = 'Differences allowed in the pre.cluster step. Should be 1 every 100 bases.If the data are only Pacbio, it is ignored'
@params['minLen'] = '145'
@params['minLen', 'description'] = 'Sequences shorter than this long are removed.'
@params['maxLen'] = '330'
@params['maxLen', 'description'] = 'Sequences longer than this are removed.'
@params['mail'] = ""
@inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
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
def next_dataset
    {'Name'=>@dataset['Name'],
     'CountTable [File]'=>File.join(@result_dir, "#{@dataset['Name']}.count.txt"),
     'PreClusteredFastaFile [File]'=>File.join(@result_dir, "#{@dataset['Name']}.preclustered.txt"),
     'TaxonomyFile [File]'=>File.join(@result_dir, "#{@dataset['Name']}.taxonomy.txt"),
     }
end
def commands
run_RApp("EzAppMothurDataCleanBatch")
end
end

if __FILE__ == $0
end
