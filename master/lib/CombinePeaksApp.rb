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
    @required_columns = ['Name','BED', 'BAM', 'refBuild', 'refFeatureFile', 'Species', 'paired','Read Count']
    @required_params = ['skipExtraChr', 'minSamples']
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '20'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '100'
    @params['scratch', "context"] = "slurm"
    @params['skipExtraChr'] = true
    @params['minSamples'] = 2
    @params['cmdOptions'] = ''
    @params['specialOptions'] = ''
    @params['name'] = 'peakCountResult'
    @params['mail'] = ""
    @modules = ["Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  
  def next_dataset
    resultPath = File.join(@result_dir, "#{@params['name']}")
    dataset = {'Name'=>@dataset['Name'],
     'Count [File,Link]'=>File.join(resultPath, "#{@dataset['Name']}_peak_counts.txt"),
     'Species'=>@dataset['Species']
     'refBuild'=>@dataset['refBuild'],
     'featureLevel'=>'peak',
     'refFeatureFile'=>@dataset['refFeatureFile'],
     'paired'=>@dataset['paired'],
     'Read Count'=>@dataset['Read Count'],
      'PeakCountResult [File]'=>resultPath
    }.merge(extract_columns(@inherit_tags))
    dataset
  end
  def commands
    run_RApp("EzAppPeakCombiner")
  end
end

if __FILE__ == $0
  
end
