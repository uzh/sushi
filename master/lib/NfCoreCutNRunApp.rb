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
    @params['cores', "context"] = "slurm"
    @params['ram'] = '100'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '200'
    @params['scratch', "context"] = "slurm"
    @params['paired'] = true
    @params['paired', "context"] = "NfCoreCutNRun"
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "referfence genome assembly"
    @params['refFeatureFile'] = 'genes.gtf'
    @params['refFeatureFile', "context"] = "NfCoreCutNRun"
    @params['peakStyle'] = ['broad', 'narrow']
    @params['peakCaller'] = ['macs2', 'seacr']
    @params['spikeinGenome'] = ['K12-MG1655', 'R64-1-1', 'BDGP6']
    @params['spikeinGenome', 'description'] = 'name of the iGenome reference for the spike-in genome. Defaulting to E. coli K12, for yeast set to R64-1-1, for fruit fly BDGP6'
    @params['normalization'] = ['Spikein', 'RPKM', 'CPM', 'BPM', 'None']
    @params['grouping'] = ''
    @params['grouping', 'description'] = 'grouping information needs to be filled in the input dataset'
    @params['controlColumn'] = ''
    @params['controlColumn', 'description'] = 'indicate the column with the control samples'
    @params['keepBams'] = false
    @params['keepBams', 'description'] = 'delete BAM and BAI files from alignment folder'
    @params['name'] = 'NfCoreCutAndRun'
    @params['pipelineVersion'] = '3.2.2'
    @params['pipelineVersion', 'description'] = 'specify pipeline version of nf-core pipeline'
    @params['cmdOptions'] = ""
    @params['mail'] = ""
    @modules = ["Dev/jdk", "Tools/BEDTools"]
  end
 def set_default_parameters
    if @params['paired']
      @required_columns<<  'Read2'
    end
  end
  def next_dataset
     report_file = File.join(@result_dir, "#{@params['name']}_results")
     multiqc_report_link = File.join(@result_dir, "#{@params['name']}_results", "04_reporting", "multiqc", "multiqc_report.html")
     cutandrun_report_link = File.join(@result_dir, "#{@params['name']}_results", "04_reporting", "00index.html")
     igv_link = "https://fgcz-gstore.uzh.ch/projects/#{report_file}/igv_session.html"
     {'Name'=>@params['name'],
      'CutAndRun_Result [File]'=>report_file,
      'Multiqc [Link]'=>multiqc_report_link,
      'CutAndRun_Report [Link]'=>cutandrun_report_link,
      'IGV [Link]'=>igv_link,
      'Species'=>(dataset = @dataset.first and dataset['Species']),
      'refBuild'=>@params['refBuild']
    }
  end
  def commands
    run_RApp('EzAppNfCoreCutAndRun')
  end
end
