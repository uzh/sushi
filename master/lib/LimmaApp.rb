#!/usr/bin/env ruby
# encoding: utf-8
Version = '20190131-041245'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class LimmaApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Limma'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Differential_Expression'
    @description =<<-EOS
    Empirical analysis of digital gene expression data in R<br/>
<a href='https://bioconductor.org/packages/release/bioc/html/limma.html'>limma</a><br/>
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
    @params['modelMethod'] = ['limma-trend', 'voom']
    @params['modelMethod', 'description'] = 'Which model method to use'
    @params['grouping'] = '' ### TODO: this should be filled by a column selector that allows to select a column with the tag 'Factor'
    @params['sampleGroup'] = '' ## TODO: this should be a value from the selected column
    @params['sampleGroup', 'description'] = 'sampleGroup should be different from refGroup'
    @params['refGroup'] = '' ## TODO: this should be a value from the selected column
    @params['refGroup', 'description'] = 'refGroup should be different from sampleGroup'
    @params['onlyCompGroupsHeatmap'] = ['true', 'false']
    @params['onlyCompGroupsHeatmap', 'description'] = 'Only show the samples from comparison groups in heatmap'
    @params['runGO'] = ['true', 'false']
    @params['grouping2'] = ''
    @params['grouping2', 'description'] =  'specify the column name of your secondary co-variate  (factor or numeric, 
    assuming there is one). Ensure the 
    column name is in the format "NAME [Factor]" or "NAME [Numeric]"'
    @params['priorCount'] = 10
    @params['priorCount', "description"] = "counts to be added to shrink estimated log2 ratios"
    @params['transcriptTypes'] = ''
    @params['transcriptTypes', 'multi_selection'] = true
    @params['transcriptTypes', 'selected'] = 0
    @params['specialOptions'] = ''
    @params['expressionName'] = ''
    @params['mail'] = ""
    @params['Rversion'] = ["Dev/R/4.4.0", "Dev/R/4.3.2", "Dev/R/4.3.0", "Dev/R/4.2.2", "Dev/R/4.2.0", "Dev/R/4.1.2", "Dev/R/4.0.4", "Dev/R/4.0.3", "Dev/R/4.0.1"]
    @inherit_columns = ["Order Id"]
  end
   def preprocess
    @random_string = (1..12).map{[*('a'..'z')].sample}.join
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
    command = "module load #{@params["Rversion"]}\n"
    command << run_RApp("EzAppLimma")
  end
end

if __FILE__ == $0

end
