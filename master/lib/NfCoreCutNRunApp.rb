#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class NfCoreCutNRunApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'NfCoreCutNRun'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'CutAndRun'
    @description =<<-EOS
    A CUT&RUN processing pipeline from NF-Core. <br/>
   <a href='https://nf-co.re/cutandrun'>NF-Core CUT&RUN</a>
EOS
    @required_columns = ['Name','Read1','Read2','Species']
    @required_params = ['refBuild', 'peakStyle', 'grouping', 'controlColumn', 'normalization', 'spikeinGenome']
    @params['cores'] = '8'
    @params['ram'] = '100'
    @params['scratch'] = '200'
    @params['paired'] = true
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['peakStyle'] = ['broad', 'narrow']
    @params['peakCaller'] = ['macs2', 'seacr']
    @params['spikeinGenome'] = ['K12-MG1655', 'R64-1-1', 'BDGP6']
    @params['spikeinGenome', 'description'] = 'name of the iGenome reference for the spike-in genome. Defaulting to E. coli K12, for yeast set to R64-1-1, for fruit fly BDGP6'
    @params['normalization'] = ['Spikein', 'RPKM', 'CPM', 'BPM', 'None']
    @params['grouping'] = ''
    @params['grouping', 'description'] = 'grouping information needs to be filled in the input dataset'
    @params['controlColumn'] = ''
    @params['controlColumn', 'description'] = 'indicate the column with the control samples'
    @params['name'] = 'NfCoreCutAndRun'
    @params['cmdOptions'] = ""
    @params['mail'] = ""
    @modules = ["Dev/jdk"]
  end
 def set_default_parameters
    if @params['paired']
      @required_columns<<  'Read2'
    end
  end
  def next_dataset
     report_file = File.join(@result_dir, @params['name'])
     report_link = File.join(report_file, 'multiqc_report.html')
     {'Name'=>@params['name'],
      'Result [File]'=>report_file,
      'Report [Link]'=>report_link,
      'Species'=>(dataset = @dataset.first and dataset['Species']),
      'refBuild'=>@params['refBuild']
    }
  end
  def commands
    run_RApp('EzAppNfCoreCutAndRun')
  end
end
