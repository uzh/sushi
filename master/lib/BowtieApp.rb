#!/usr/bin/env ruby
# encoding: utf-8
Version = '20160215-003903'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class BowtieApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Bowtie'
    @analysis_category = 'Map'
    @description =<<-EOS
Fast and memory-efficient short read aligner<br/>
<a href='http://bowtie-bio.sourceforge.net/index.shtml'>manual</a>
    EOS
    @required_columns = ['Name','Read1','Species']
    @required_params = ['refBuild','paired']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '16'
    @params['scratch'] = '100'
    @params['refBuild'] = ref_selector
    @params['paired'] = false
    @params['cmdOptions'] = ''
    @params['trimAdapter'] = true
    @params['trimLeft'] = 0
    @params['trimRight'] = 0
    @params['minTailQuality'] = 0
    @params['minAvgQuality'] = 10
    @params['minReadLength'] = 20
    @params['specialOptions'] = ''
    @params['mail'] = ""
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
     'BAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam"),
     'BAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam.bai"),
     'IGV Starter [Link]'=>File.join(@result_dir, "#{@dataset['Name']}-igv.jnlp"),
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'paired'=>@params['paired'],
     'Read Count'=>@dataset['Read Count'],
     'IGV Starter [File]'=>File.join(@result_dir, "#{@dataset['Name']}-igv.jnlp"),
     'IGV Session [File]'=>File.join(@result_dir, "#{@dataset['Name']}-igv.xml"),
     'PreprocessingLog [File]'=>File.join(@result_dir, "#{@dataset['Name']}_preprocessing.log")      
    }.merge(extract_column("Factor")).merge(extract_column("B-Fabric"))
  end
  def commands
    run_RApp("EzAppBowtie")
  end
end

if __FILE__ == $0
  #run BowtieApp

end

