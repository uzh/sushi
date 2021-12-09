#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class CountQCApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'CountQC'
    @analysis_category = 'QC'
    @description =<<-EOS
    Run test module in the tool Model-based Analysis of Genome-wide CRISPR-Cas9 Knockout (<a href='https://sourceforge.net/p/mageck/wiki/Home/'>MAGeCK</a>)
    EOS
    @params['process_mode'] = 'DATASET'
    @required_columns = ['Name','Count']
    @required_params = []
    # optional params
    @params['cores'] = ['1']
    @params['ram'] = ['4']
    @params['scratch'] = ['10']
    @params['name'] = 'MAGeCK_Test'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Dev/R"]
    @inherit_tags = ["B-Fabric"]
  end
  def preprocess
  end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@params['name'],
     'Species'=>(dataset = @dataset.first and dataset['Species']),
     'refBuild'=>@params['refBuild'],
     'Static Report [Link]'=>report_link,
     'Live Report [Link]'=>"#{SHINY_EXPLORE_COUNTS}?data=#{report_file}/counts-#{@random_string}-EzResult.RData",
     'Report [File]'=>report_file,
    }.merge(extract_columns(@inherit_tags))
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('refFeatureFile')
      @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    end
    @params['transcriptTypes'] = @dataset[0]['transcriptTypes'].to_s.split(',')
  end
  def commands
    run_RApp("EzAppCountQC")
  end
end

if __FILE__ == $0

end

