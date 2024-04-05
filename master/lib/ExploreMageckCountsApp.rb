#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class ExploreMageckCountsApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'ExploreMageckCounts'
    @analysis_category = 'QC'
    @description =<<-EOS
Quality control after counting sgRNAs  with Mageck<br/>
    EOS
    @params['process_mode'] = 'DATASET'
    @required_columns = ['Name','Count', 'Species']
    @required_params = []
    # optional params
    @params['cores'] = ['1', '2']
    @params['ram'] = ['4', '7', '50']
    @params['scratch'] = ['10', '20']
    @params['name'] = 'ExploreMageckCounts'
    #@params['refBuild'] = ref_selector
    #@params['refFeatureFile'] = 'genes.gtf'
    #@params['featureLevel'] = ['gene', 'isoform']
    @params['normMethod'] = 'logMean'
    #@params['runGO'] = ['true', 'false']
    @params['backgroundExpression'] = 5
    @params['backgroundExpression', "description"] = "counts to be added to shrink estimated log2 ratios"
    #@params['transcriptTypes'] = ''
    #@params['transcriptTypes', 'multi_selection'] = true
    #@params['transcriptTypes', 'selected'] = 0
    @params['specialOptions'] = ''
    @params['expressionName'] = 'sample1'
    @params['mail'] = ""
    @modules = ["Dev/R"]
    @inherit_tags = ["B-Fabric"]
  end
  def preprocess
    @random_string = (1..12).map{[*('a'..'z')].sample}.join
  end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@params['name'],
     'Species'=>(dataset = @dataset.first and dataset['Species']),
     #'refBuild'=>@params['refBuild'],
     'libName'=>(dataset = @dataset.first and dataset['libName']),
     'Static Report [Link]'=>report_link,
     'Live Report [Link]'=>"#{SHINY_EXPLORE_COUNTS}?data=#{report_file}/counts-#{@random_string}-EzResult.RData",
     'Report [File]'=>report_file,
    }.merge(extract_columns(@inherit_tags))
  end
  def set_default_parameters
    #@params['refBuild'] = @dataset[0]['refBuild']
    #if dataset_has_column?('refFeatureFile')
    #  @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    #end
    #@params['transcriptTypes'] = @dataset[0]['transcriptTypes'].to_s.split(',')
  end
  def commands
    run_RApp("EzAppExploreMageckCounts")
  end
end

if __FILE__ == $0

end

