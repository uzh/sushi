#!/usr/bin/env ruby
# encoding: utf-8
Version = '20180412-100213'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class DESeq2App < SushiFabric::SushiApp
  def initialize
    super
    @name = 'DESeq2'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Differential_Expression'
    @description =<<-EOS
    Differential gene expression analysis based on the negative binomial distribution<br/>
<a href='https://bioconductor.org/packages/release/bioc/html/DESeq2.html'>DESeq2</a><br/>
    EOS
    @required_columns = ['Name','Count', 'Species', 'refBuild', 'featureLevel', 'refFeatureFile']
    @required_params = ['grouping', 'sampleGroup', 'refGroup']
    # optional params
    @params['cores'] = '1'
    @params['ram'] = '2'
    @params['scratch'] = '10'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['featureLevel'] = ['gene', 'isoform']
    @params['grouping'] = ''
    @params['sampleGroup'] = ''
    @params['sampleGroup', 'description'] = 'sampleGroup should be different from refGroup'
    @params['refGroup'] = ''
    @params['refGroup', 'description'] = 'refGroup should be different from sampleGroup'
    #@params['normMethod'] = 'logMean'
    @params['runGO'] = ['true', 'false']
    @params['grouping2'] = ''
    @params['grouping2', 'description'] = 'specify the column name of your secondary factor --only in case your experiment has a second factor, that should enter the linear model for differential expression!'
    @params['backgroundExpression'] = 10
    @params['backgroundExpression', "description"] = "counts to be added to shrink estimated log2 ratios"
    @params['transcriptTypes'] = ''
    @params['transcriptTypes', 'multi_selection'] = true
    @params['transcriptTypes', 'selected'] = 0
    @params['specialOptions'] = ''
    @params['expressionName'] = ''
    @params['mail'] = ""
    @modules = ["Tools/samtools", "Tools/GFOLD", "Dev/PhantomJS", "Dev/R", "Tools/sambamba"]
  end
  def preprocess
    @random_string = (1..12).map{[*('a'..'z')].sample}.join
    @name = "#{@name}_#{@params['sampleGroup']}--over--#{@params['refGroup']}"
  end
  def next_dataset
    @comparison = "#{@params['sampleGroup']}--over--#{@params['refGroup']}"
    @params['comparison'] = @comparison
    @params['name'] = @comparison
    report_file = File.join(@result_dir, "#{@params['comparison']}")
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@comparison,
     'Species'=>(dataset = @dataset.first and dataset['Species']),
     'refBuild'=>@params['refBuild'],
     'Static Report [Link]'=>report_link,
     'Live Report [Link]'=>"#{SHINY_EXPLORE_DE}?data=#{report_file}/result-#{@comparison}-#{@random_string}-EzResult.RData",
     'Report [File]'=>report_file,
    }
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('refFeatureFile')
      @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    end
    @params['transcriptTypes'] = @dataset[0]['transcriptTypes'].to_s.split(',')
  end
  def commands
    run_RApp("EzAppDeseq2")
  end
end

if __FILE__ == $0

end
