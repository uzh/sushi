#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class NfCoreAtacSeqApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'NfCoreAtacSeq'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'ATAC'
    @description =<<-EOS
    A ATAC-seq processing pipeline from NF-Core. <br/>
   <a href='https://nf-co.re/atacseq'>NF-Core ATAC-seq</a>
EOS
    @required_columns = ['Name','Read1','Read2','Species']
    @required_params = ['refBuild', 'peakStyle', 'grouping']
    @params['cores'] = '8'
    @params['ram'] = '100'
    @params['scratch'] = '200'
    @params['paired'] = true
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['peakStyle'] = ['broad', 'narrow']
    @params['varStabilizationMethod'] = ['vst', 'rlogTransf']
    @params['grouping'] = ''
    @params['grouping', 'description'] = 'grouping information needs to be filled in the input dataset'
    @params['runTwoGroupAnalysis'] = false
    @params['runTwoGroupAnalysis', 'description'] = 'perform all two group analysis based on the grouping information'
    @params['keepBams'] = false
    @params['keepBams', 'description'] = 'delete BAM and BAI files from the result folder'
    @params['name'] = 'NfCoreAtacSeq'
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
     ## the line below demonstrates that access to @dataset does not work as expected when in dataset_mode
     #foo = @dataset['Name']
     report_file = File.join(@result_dir, "#{@params['name']}_result")
     multiqc_link = File.join(@result_dir, "#{@params['name']}_result", "multiqc", "#{@params['peakStyle']}_peak", "multiqc_report.html")
     ataqv_link = File.join(@result_dir, "#{@params['name']}_result", "bwa/merged_library/ataqv", "#{@params['peakStyle']}_peak", "html/index.html")
     igv_link = File.join(@result_dir, "#{@params['name']}_result", "igv", "#{@params['peakStyle']}_peak", "igv_session.html")
     dataset = {'Name'=>@params['name'],
     'Species'=>(dataset = @dataset.first and dataset['Species']),
     'refBuild'=>@params['refBuild'],
     'ATAC_Result [File]'=>report_file,
     'Multiqc [Link]'=>multiqc_link,
     'Ataqv [Link]'=> ataqv_link,
     'IGV [Link]'=>igv_link
     }
    if @params['runTwoGroupAnalysis']
      diffReport_link = File.join(@result_dir, "#{@params['name']}_results", "diffpeak_analysis", "DifferentialPeaks.html")
      dataset['DifferentialPeaks [Link]'] = diffReport_link
    end
    dataset
  end
  def grandchild_datasets
    ## TODO ... this is a workaround to get grandchild_data that is per sample; I need the grandchild data per sample in _one_ dataset
    my_dataset = []
    grandchild_data = []
    @dataset_hash.each_with_index do |row, i|
      my_dataset = Hash[*row.map{|key,value| [key.gsub(/\[.+\]/,'').strip, value]}.flatten]
      grandchild_data << {
        'Name'=>my_dataset['Name'],
        
        'Count [File]'=>File.join(@result_dir, "#{@params['name']}_result", "bwa/merged_library/macs2/broad_peak/consensus", "#{my_dataset['Name']}.txt"),
        'BigWig [File]'=>File.join(@result_dir, "#{@params['name']}_result", "bwa/merged_library/bigwig", "#{my_dataset['Name']}.bigWig"),
        'Species'=>my_dataset['Species'],
        'refBuild'=>@params['refBuild'],
        'featureLevel'=>@params['peakStyle']
      } ## TODO I would like to merge from the parent .merge(extract_columns(@inherit_tags))
    end
    grandchild_data
  end
  def commands
    run_RApp('EzAppNfCoreAtacSeq')
  end
end
