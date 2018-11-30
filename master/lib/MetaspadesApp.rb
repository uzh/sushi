#!/usr/bin/env ruby
# encoding: utf-8
Version = '20181127-155223'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class MetaspadesApp < SushiFabric::SushiApp
def initialize
super
@name = 'Metaspades'
@analysis_category = 'Metagenomics'
@description =<<-EOS
Denovo metagenomics assembly with Metaspades, gene prediction  with Prodigal and annotation with InterProScan (soon also Diamond). 
<a href='http://cab.spbu.ru/software/spades/'>http://cab.spbu.ru/software/spades/</a>
<a href='https://github.com/hyattpd/prodigal/wiki'>https://github.com/hyattpd/prodigal/wiki</a>
<a href='https://github.com/bbuchfink/diamond'>https://github.com/bbuchfink/diamond</a>
<a href='https://github.com/ebi-pf-team/interproscan/wiki/InterProScan5RC4'>https://github.com/ebi-pf-team/interproscan/wiki/InterProScan5RC4</a>
  EOS
@required_columns = ['Name', 'Read1']
@required_params = ['metaspadesKmerList','diamondEvalue','diamondMaxSeqs']
@params['cores'] = '1'
@params['ram'] = '8'
@params['scratch'] = '10'
@params['metaspadesKmerList'] = '69,79,89'
@params['metaspadesKmerList', 'description'] = 'Comma-separated list of k-mer for the assembly.'
@params['diamondEvalue'] = '0.05'
@params['diamondEvalue', 'description'] = 'Blast e-value cut-off.'
@params['diamondMaxSeqs'] = '30'
@params['diamondMaxSeqs', 'description'] = 'Blast maximum number of sequences to report.'
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
     'contigFile [File]' => File.join(@result_dir, "#{@dataset['Name']}.contigs.fasta"),
     'prodigalPredictionFile [File]' => File.join(@result_dir, "#{@dataset['Name']}.prodigalAnnotation.gff"),
     'interproscanFile [File]' => File.join(@result_dir, "#{@dataset['Name']}.annotatedProteins.gff"),
     'krakenLabelsFile [File]' => File.join(@result_dir, "#{@dataset['Name']}.kraken.labels")
}
end
def commands
run_RApp("EzAppMetaspades")
end
end

if __FILE__ == $0
end
