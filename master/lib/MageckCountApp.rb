#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class MageckCountApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'MageckCount'
    @params['process_mode'] = 'SAMPLE'
    @analysis_category = 'Misc'
    @description =<<-EOS
    
Model-based Analysis of Genome-wide CRISPR-Cas9 Knockout (<a href='https://sourceforge.net/p/mageck/wiki/Home/'>MAGeCK</a>) is a computational tool to identify important genes from the recent genome-scale CRISPR-Cas9 knockout screens
(or GeCKO) technology. MAGeCK is developed by Wei Li and Han Xu from Dr. Xiaole Shirley Liu's lab at Dana-Farber Cancer Institute, 
and is being actively updated by Wei Li lab from Children's National Medical Center. 
This app counts the specified sgRNAs based on raw data fastq files.
EOS

    @required_columns = ['Name','Read1']
    @required_params = ['libName']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '100'

    @params['name'] = 'MAGeCK_CountResult'
    @params['libName'] = ''
    @params['libName'] = {'select'=>''}
    Dir["/srv/GT/databases/GEML/sgRNA_Libs/*"].sort.select{|lib| File.directory?(lib)}.each do |dir|
      @params['dictPath'][File.basename(dir)] = File.basename(dir)
    end
    ## additional commands
    @params['mail'] = ""
    @modules = ["Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def preprocess
    #if @params['paired']
    #  @required_columns << 'Read2'
    #end
  end
 def set_default_parameters
    #@params['paired'] = dataset_has_column?('Read2')
  end
  def next_dataset
    report_file = File.join(@result_dir,"#{@dataset['Name']}")
    count_file = File.join(report_file, "#{@dataset['Name']}.count.txt")
    log_file = File.join(report_file, "#{@dataset['Name']}.log")
    {'Name'=>@dataset['Name'],
     'Count [File]'=>count_file,
     'Log [File]'=>log_file,
     'Read Count'=>@dataset['Read Count'] 
    }
  end
  def commands
    run_RApp("EzAppMageckCount")
  end
end

if __FILE__ == $0
  run MageckCountApp
end
