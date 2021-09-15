#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class InlineBarcodeDmxApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'InlineBarcodeDmx'
    @params['process_mode'] = 'SAMPLE'
    @analysis_category = 'Prep'
    @description =<<-EOS
    Simple R based app which performs demultiplexing of inline barcodes, read trimming and read filtering based on remaining read length
EOS
    @required_columns = ['Name','Read1','indexFile']
    @required_params = ['name']
    @params['cores'] = '1'
    @params['ram'] = '30'
    @params['scratch'] = '200'
    @params['name'] = 'InLineDmx_Result'
    @params['barcodePos'] = '8'
    @params['rightPattern'] = 'GTGTCAGTCACTTCCAG'
    @params['maxMismatch'] = '1'
    @params['minReadLength'] = '20'
    @params['maxReadLength'] = '50'
    @params['leftLinkerSize'] = '5'

    @params['cmdOptions'] = ""
    @params['mail'] = ""
    @modules = ["Dev/R"]
  end
# def set_default_parameters
#    @params['paired'] = dataset_has_column?('Read2')
#  end
#  def preprocess
#    if @params['paired']
#      @required_columns<<  'Read2'
#    end
#  end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])    
    {'Name'=>@params['name'],
     'Report [File]'=>report_file,
    }
  end
  def commands
    run_RApp("EzAppInlineBarcodeDmx")
  end
end

if __FILE__ == $0
  
end
