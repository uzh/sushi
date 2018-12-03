#!/usr/bin/env ruby
# encoding: utf-8
Version = '20180905-113400'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class UparseApp < SushiFabric::SushiApp
def initialize
super
@name = 'Uparse'
@analysis_category = 'Metagenomics'
@description =<<-EOS
OTU-based metagenomics analysis with Uparse. 
Please make sure that the input files are from the same technology.
<a href='https://drive5.com/uparse/'>https://drive5.com/uparse/</a>
  EOS
@params['process_mode'] = 'DATASET'
@required_columns = ['Name', 'Read1']
@required_params = ['fastqErrorMax']
@params['cores'] = '1'
@params['ram'] = '8'
@params['scratch'] = '10'
@params['Group'] = ['false','true']
@params['Group', 'description'] = 'Is there a design matrix for the experiment? '
@params['fastqErrorMax'] = '1'
@params['fastqErrorMax', 'description'] = 'Max fastq error rate (https://drive5.com/usearch/manual/cmd_fastq_filter.html).'
@params['mail'] = ""
@params['Name'] = "Uparse"
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
     nds = {'Name'=>@params['name']}
     nds['OTUsToTaxonomyFile [File]'] = File.join(@result_dir, "#{@params['Name']}.OTUs.to.tax.txt")
     nds['OTUsCountTable [File]'] = File.join(@result_dir, "#{@params['Name']}.OTUs.count.txt")
      if @params['Group']
     nds['sampleDescriptionFile [File]'] = File.join(@result_dir, "#{@params['Name']}.designMatrix.txt")
      end
     nds
end
def commands
run_RApp("EzAppUparse")
end
end

if __FILE__ == $0
end
