#!/usr/bin/env ruby
# encoding: utf-8
Version = '20180905-113400'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class DADA2Step2DatasetApp < SushiFabric::SushiApp
def initialize
super
@name = 'DADA2Step2Dataset'
@analysis_category = 'Metagenomics'
@description =<<-EOS
Data preprocssing with DADA2. Please make sure that the input files are from the same technology and adjust minLen and maxLen accordingly.
<a href='https://DADA2.org/wiki/MiSeq_SOP'>https://DADA2.org/wiki</a>
  EOS
@params['process_mode'] = 'DATASET'
@required_columns = ['Name', 'RObjectWithSeqTab']
@required_params = ['database']
@params['cores'] = '2'
@params['ram'] = '8'
@params['scratch'] = '10'
@params['database'] = ['silva','RDP','greenGenes']
@params['database', 'description'] = '16S database to use for taxonomic assignment.'
@params['mail'] = ""
@inherit_tags = ['B-Fabric', 'Characteristic', 'Mock','Group']
@modules = ['Dev/R']
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
     nds = {'Name'=>@params['Name']}
      if @params['Group']
     nds['sampleDescriptionFile [File]'] = File.join(@result_dir, "#{@params['Name']}.designMatrix.txt")
      end
     nds['OTUsToTaxonomyFile [File]'] = File.join(@result_dir, "#{@params['Name']}.OTUsToTax.txt")
     nds['OTUsCountTable [File]'] = File.join(@result_dir, "#{@params['Name']}.OTUsToCount.txt")
     nds
end
def commands
run_RApp("EzAppDADA2Step2Dataset")
end
end

if __FILE__ == $0
end
