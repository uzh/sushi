#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class AtacSeqApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'AtacSeq'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'GeneRegulation'
    @description =<<-EOS
    A ATAC-seq processing pipeline from NF-Core. <br/>
   <a href='https://nf-co.re/atacseq'>NF-Core ATAC-seq</a>
EOS
    @required_columns = ['Name','Read1', 'Species']
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
    @required_params = ['refBuild', 'peakStyle']
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '100'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '200'
    @params['scratch', "context"] = "slurm"
    @params['paired'] = true
    @params['paired', "context"] = "AtacSeq"
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "reference genome assembly"
    @params['refFeatureFile'] = 'genes.gtf'
    @params['refFeatureFile', "context"] = "AtacSeq"
    @params['peakStyle'] = ['broad', 'narrow']
    @params['varStabilizationMethod'] = ['vst', 'rlogTransf']
    @params['grouping'] = 'Condition'
    @params['grouping', 'description'] = 'dataset column with the sample groups (replicates are merged by nf-core); samples without a group are processed as their own group'
    @params['keepBams'] = false
    @params['keepBams', 'description'] = 'keep the BAM and BAI files in the result folder'
    @params['name'] = 'AtacSeq'
    @params['pipelineVersion'] = '2.1.2'
    @params['pipelineVersion', 'description'] = 'specify pipeline version of nf-core pipeline'
    @params['qcMode'] = false
    @params['qcMode', 'description'] = 'QC run: use only the first qcReadsPerSample reads (pairs) of each sample'
    @params['qcReadsPerSample'] = ['10000000', '5000000', '20000000', '50000000']
    @params['qcReadsPerSample', 'description'] = 'number of reads (pairs) per sample in QC mode; samples with fewer reads are used completely'
    @params['cmdOptions'] = ""
    @params['cmdOptions', "context"] = "AtacSeq"
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
