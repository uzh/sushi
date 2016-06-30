#!/usr/bin/env ruby
# encoding: utf-8
Version = '20160215-003128'

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
    @params['refGroup'] = ''
    #@params['normMethod'] = 'logMean'
    @params['runGO'] = ['false', 'true']
    @params['expressionName'] = ''
    @params['specialOptions'] = ''
    @params['mail'] = ""
  end
  def preprocess
    @random_string = (1..12).map{[*('a'..'z')].sample}.join
  end
  def next_dataset
    @comparison = "#{@params['sampleGroup']}--over--#{@params['refGroup']}"
    @params['comparison'] = @comparison
    @params['name'] = @comparison
    report_file = File.join(@result_dir, "#{@name}--#{@params['name']}")
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@comparison,
     'Species'=>@dataset['Species'],
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
  end
  def commands
    run_RApp("EzAppDeseq2")
  end
end

if __FILE__ == $0

end

