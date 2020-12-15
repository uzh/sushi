#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-095333'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class ScaterApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'ScaterApp'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
    Maps all read files specified by the dataset file generates stats and expression counts with featureCounts
EOS
    @required_columns = ['Name','CountDataset', 'CountMatrix', 'Species']
    @required_params = ['refBuild']
    # optional params
    @params['cores'] = '1'
    @params['ram'] = '7'
    @params['scratch'] = '10'
    @params['name'] = 'Scater_QC'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['featureLevel'] = ['gene', 'isoform']
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Dev/R"]
  end
  def preprocess
    @random_string = (1..12).map{[*('a'..'z')].sample}.join
  end
  def next_dataset
    report_file = File.join(@result_dir, "#{@dataset['Name']}_ScaterQC")
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@dataset['Name'],
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'Static Report [Link]'=>report_link,
     'Live Report [Link]'=>"#{SHINY_SCATER}?data=#{report_file}/counts-#{@random_string}-EzResult.RData",
     'Report [File]'=>report_file,
    }
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('refFeatureFile')
      @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    end
  end
  def commands
    run_RApp("EzAppScater")
  end
end

