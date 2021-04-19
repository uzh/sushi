#!/usr/bin/env ruby
# encoding: utf-8
Version = '20201127-101634'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class EdgeRApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'EdgeR'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Differential_Expression'
    @description =<<-EOS
    Empirical analysis of digital gene expression data in R<br/>
<a href='https://bioconductor.org/packages/release/bioc/html/edgeR.html'>edgeR</a><br/>
    EOS
    @required_columns = ['Name','Count', 'Species', 'refBuild', 'featureLevel', 'refFeatureFile']
    @required_params = ['grouping', 'sampleGroup', 'refGroup']
    # optional params
    @params['cores'] = ['4', '2']
    @params['ram'] = ['16', '12', '8']
    @params['scratch'] = ['10', '20']
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['featureLevel'] = ['gene', 'isoform']
    @params['testMethod'] = ['glm', 'exactTest']
    @params['deTest'] = ['QL', 'LR']
    @params['deTest', 'description'] = 'This option only works for glm method. Quasi-likelihood (QL) F-test or likelihood ratio (LR) test. LR is prefered for single-cell data.'
    @params['grouping'] = '' ## Note: this will be a selector defined by Factor tagged column
    @params['sampleGroup'] = '' ## Note: this will be a selector defined by Factor tagged column
    @params['sampleGroup', 'description'] = 'sampleGroup should be different from refGroup'
    @params['sampleGroupBaseline'] = '' ### Note: this will be a selector defined by Factor tagged column
    @params['sampleGroupBaseline', 'description'] = 'select the baseline for sampleGroup if you have'
    @params['refGroup'] = '' ## Note: this will be a selector defined by Factor tagged column
    @params['refGroup', 'description'] = 'refGroup should be different from sampleGroup'
    @params['refGroupBaseline'] = '' ## Note: this will be a selector defined by Factor tagged column
    @params['refGroupBaseline', 'description'] = 'select the baseline for refGroup if you have'
    @params['onlyCompGroupsHeatmap'] = ['false', 'true']
    @params['onlyCompGroupsHeatmap', 'description'] = 'Only show the samples from comparison groups in heatmap'
    @params['normMethod'] = ['TMM', 'RLE', 'upperquartile', 'none']
    @params['normMethod', 'description'] = "see http://bioconductor.org/packages/edgeR/"
    @params['runGO'] = ['true', 'false']
    @params['grouping2'] = ''
    @params['grouping2', 'description'] =  'specify the column name of your secondary co-variate  (factor or numeric, 
    assuming there is one). Ensure the 
    column name is in the format "NAME [Factor]" or "NAME [Numeric]"'
    @params['backgroundExpression'] = 10
    @params['backgroundExpression', "description"] = "counts to be added to shrink estimated log2 ratios"
    @params['transcriptTypes'] = ''
    @params['transcriptTypes', 'multi_selection'] = true
    @params['transcriptTypes', 'selected'] = 0
    @params['specialOptions'] = ''
    @params['expressionName'] = ''
    @params['mail'] = ""
    @params['Rversion'] = ["Dev/R/4.0.4", "Dev/R/4.0.3", "Dev/R/4.0.1"]
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
    command = "module load #{@params["Rversion"]}\n"
    command << run_RApp("EzAppEdger")
  end
end

if __FILE__ == $0

end
