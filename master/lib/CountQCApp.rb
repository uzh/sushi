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
Quality control after counting reads<br/>
    EOS
    @params['process_mode'] = 'DATASET'
    @required_columns = ['Name','Count', 'Species', 'refBuild', 'featureLevel', 'refFeatureFile']
    @required_params = []
    # optional params
    @params['cores'] = ['1', '2']
    @params['cores', "context"] = "slurm"
    @params['ram'] = ['8', '12', '50']
    @params['ram', "context"] = "slurm"
    @params['scratch'] = ['10', '20']
    @params['scratch', "context"] = "slurm"
    @params['name'] = 'Count_QC'
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "referfence genome assembly"
    @params['refFeatureFile'] = 'genes.gtf'
    @params['refFeatureFile', "context"] = "CountQC"
    @params['featureLevel'] = ['gene', 'isoform']
    @params['featureLevel', "context"] = "CountQC"
    @params['normMethod'] = 'logMean'
    @params['runGO'] = ['true', 'false']
    @params['backgroundExpression'] = 10
    @params['backgroundExpression', "description"] = "counts to be added to shrink estimated log2 ratios"
    @params['topGeneSize'] = 100
    @params['selectByFtest'] = ['false', 'true']
    @params['transcriptTypes'] = ''
    @params['transcriptTypes', 'multi_selection'] = true
    @params['transcriptTypes', 'selected'] = 0
    @params['specialOptions'] = ''
    @params['expressionName'] = ''
    @params['mail'] = ""
    @modules = ["Dev/R"]
    @inherit_columns = ["Order Id"]
  end
  def preprocess
    @random_string = (1..12).map{[*('a'..'z')].sample}.join
  end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@params['name'],
     'Species'=>(dataset = @dataset.first and dataset['Species']),
     'refBuild'=>@params['refBuild'],
     'Static Report [Link]'=>report_link,
     'Live Report [Link]'=>"http://fgcz-shiny.uzh.ch/exploreDE?data=#{report_file}/counts-#{@random_string}-EzResult.RData",
     'Report [File]'=>report_file,
    }.merge(extract_columns(colnames: @inherit_columns))
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

