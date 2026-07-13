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
    @analysis_category = 'GeneRegulation'
    @description =<<-EOS
    A ATAC-seq processing pipeline from NF-Core. <br/>
   <a href='https://nf-co.re/atacseq'>NF-Core ATAC-seq</a>
EOS
    @required_columns = ['Name','Read1', 'Species']
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
    @required_params = ['refBuild', 'peakStyle', 'grouping']
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '100'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '200'
    @params['scratch', "context"] = "slurm"
    @params['paired'] = true
    @params['paired', "context"] = "NfCoreAtacSeq"
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "referfence genome assembly"
    @params['refFeatureFile'] = 'genes.gtf'
    @params['refFeatureFile', "context"] = "NfCoreAtacSeq"
    @params['peakStyle'] = ['broad', 'narrow']
    @params['varStabilizationMethod'] = ['vst', 'rlogTransf']
    @params['grouping'] = ''
    @params['grouping', 'description'] = 'grouping information needs to be filled in the input dataset'
    @params['runTwoGroupAnalysis'] = false
    @params['runTwoGroupAnalysis', 'description'] = 'perform all two group analysis based on the grouping information'
    @params['keepBams'] = false
    @params['keepBams', 'description'] = 'delete BAM and BAI files from the result folder'
    @params['name'] = 'NfCoreAtacSeq'
    @params['pipelineVersion'] = '2.1.2'
    @params['pipelineVersion', 'description'] = 'specify pipeline version of nf-core pipeline'
    @params['cmdOptions'] = ""
    @params['cmdOptions', "context"] = "NfCoreAtacSeq"
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
     ##igv_link = File.join(@result_dir, "#{@params['name']}_result", "igv", "#{@params['peakStyle']}_peak", "igv_session.html")
     igv_link = "https://fgcz-gstore.uzh.ch/projects/#{report_file}/igv_session.html"
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
    grandchild_dataset = []
    rows = @dataset.is_a?(Array) ? @dataset : []
    return grandchild_dataset if rows.empty?
    @params['grandchildName'] = "details" ## TODO: order name should be kept
    rows.each_with_index do |row, i|
      sample = Hash[*row.map{|key,value| [key.gsub(/\[.+\]/,'').strip, value]}.flatten]
      grandchild_dataset << {
        'Name'=>sample['Name'],
        'Count [Link]'=>File.join(@result_dir, "#{@params['name']}_result", "bwa/merged_library/macs2", "#{@params['peakStyle']}_peak", "consensus", "#{sample['Name']}.txt"),
        'BigWig' =>File.join(@result_dir, "#{@params['name']}_result", "bwa/merged_library/bigwig", "#{sample['Name']}.bigWig"),
        'Species'=>sample['Species'],
        'refBuild'=>@params['refBuild'],
        'featureLevel'=>"peaks"
      }.merge(extract_columns(tags: @inherit_tags, sample_name: sample['Name']))
    end
    grandchild_dataset
  end
  def commands
    run_RApp('EzAppNfCoreAtacSeq')
  end
end
