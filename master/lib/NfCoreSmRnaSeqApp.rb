#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class NfCoreSmRnaSeqApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'NfCoreSmRnaSeq'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Count'
    @description =<<-EOS
    A smallRNA Seq processing pipeline from NF-Core. <br/>
   <a href='https://nf-co.re/smrnaseq'>NF-Core smallRNA-seq</a>
EOS
    @required_columns = ['Name','Read1', 'Species']
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
    @required_params = ['referenceGenome','mirtraceSpecies','pipelineVersion']
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '200'
    @params['referenceGenome'] = 'hg38'
    @params['referenceGenome', 'description'] = 'all supported genomes are listed under https://github.com/nf-core/rnaseq/blob/e049f51f0214b2aef7624b9dd496a404a7c34d14/conf/igenomes.config'
    @params['mirtraceSpecies'] = 'hsa'
    @params['mirtraceSpecies', 'description'] = 'It should point to the 3-letter species name used by miRBase'
    @params['name'] = 'NfCoreSmRnaSeq'
    @params['pipelineVersion'] = '2.4.0'
    @params['pipelineVersion', 'description'] = 'specify pipeline version of nf-core pipeline'
    @params['cmdOptions'] = ""
    @params['mail'] = ""
    @modules = ["Dev/jdk"]
  end
  def next_dataset
     ## the line below demonstrates that access to @dataset does not work as expected when in dataset_mode
     #foo = @dataset['Name']
     report_file = File.join(@result_dir, "#{@params['name']}_result")
     multiqc_link = File.join(@result_dir, "#{@params['name']}_result", "multiqc", "#{@params['peakStyle']}_peak", "multiqc_report.html")
     dataset = {'Name'=>@params['name'],
     'Species'=>(dataset = @dataset.first and dataset['Species']),
     'referenceGenome'=>@params['referenceGenome'],
     'smRNASeq_Result [File]'=>report_file,
     'Multiqc [Link]'=>multiqc_link
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
        'Count [Link]'=>File.join(@result_dir, "#{@params['name']}_result", "mirna_quant", "#{sample['Name']}.txt"),
        'Species'=>sample['Species'],
        'refBuild'=>"", 
        'featureLevel'=>"smRNA"
      }.merge(extract_columns(tags: @inherit_tags, sample_name: sample['Name']))
    end
    grandchild_dataset
  end
  def commands
    run_RApp('EzAppNfCoreSmRnaSeq')
  end
end
