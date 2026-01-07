#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class CombinePeaksApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'PeakCombiner'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'GeneRegulation'
    @description =<<-EOS
    Determine consensus peaks of ATAC Seq or ChIP Seq data and quantify them <br/>
EOS
    @required_columns = ['Name','BED', 'BAM', 'refBuild', 'refFeatureFile', 'Species', 'paired']
    @required_params = ['skipExtraChr', 'minSamples']
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '20'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '100'
    @params['scratch', "context"] = "slurm"
    @params['refBuild'] = ref_selector
    @params['skipExtraChr'] = true
    @params['minSamples'] = 2
    @params['cmdOptions'] = ''
    @params['specialOptions'] = ''
    @params['name'] = 'peakCountResult'
    @params['mail'] = ""
    @modules = ["Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
  end  
  def next_dataset
    ds = @dataset.first
    resultPath = File.join(@result_dir, "#{@params['name']}")
    dataset = {'Name'=>@params['name'],
      'PeakCountResult [File]'=>resultPath,
      'Species'=>ds['Species'],
      'refBuild'=>@params['refBuild'],
      'IGV_Session [Link]' => File.join(resultPath, "igv_session.html")
    }.merge(extract_columns(@inherit_tags))
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
        'Count [Link]'=>File.join(@result_dir, "#{@params['name']}", "#{sample['Name']}_peak_counts.txt"),
        'Species'=>sample['Species'],
        'refBuild'=>sample['refBuild'], 
        'featureLevel'=>"peaks",
        'refFeatureFile'=>sample['refFeatureFile'],
        'BigWig'=>sample['BigWigFile'],
      }.merge(extract_columns(tags: @inherit_tags, sample_name: sample['Name']))
    end
    grandchild_dataset
  end

  
  def commands
    run_RApp("EzAppPeakCombiner")
  end
end

if __FILE__ == $0
  
end
