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
    @params['barcodePos', 'description'] = 'Start position of the inline barcode in R1'
    @params['leftLinkerSize'] = '5'
    @params['leftLinkerSize', 'description'] = 'Length of the linker next to the barcode.'
    @params['rightPattern'] = 'GTGTCAGTCACTTCCAG'
    @params['rightPattern', 'description'] = 'Sequence at the 3 prime end to get trimmed.'
    @params['maxMismatch'] = '1'
    @params['maxMismatch', 'description'] = 'Accepted number of mismatches in the read to detect the rightPattern'
    @params['minReadLength'] = '20'
    @params['minReadLength', 'description'] = 'Minimal read length after trimming. All reads below will be discarded.'
    @params['maxReadLength'] = '50'
    @params['maxReadLength', 'description'] = 'Maximal read length after trimming. All reads above will be discarded.'

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
