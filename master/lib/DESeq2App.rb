#!/usr/bin/env ruby
# encoding: utf-8
Version = '20180905-133608'

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
    @params['cores'] = ['4', '2']
    @params['ram'] = ['12', '8']
    @params['scratch'] = ['10', '20']
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['featureLevel'] = ['gene', 'isoform']
    @params['grouping'] = ''
    @params['sampleGroup'] = ''
    @params['sampleGroup', 'description'] = 'sampleGroup should be different from refGroup'
    @params['refGroup'] = ''
    @params['refGroup', 'description'] = 'refGroup should be different from sampleGroup'
    @params['onlyCompGroupsHeatmap'] = ['false', 'true']
    @params['onlyCompGroupsHeatmap', 'description'] = 'Only show the samples from comparison groups in heatmap'
    #@params['normMethod'] = 'logMean'
    @params['grouping2'] = ''
    @params['grouping2', 'description'] =  'specify the column name of your secondary co-variate (factor or numeric, 
    assuming there is one). Ensure the 
    column name is in the format "NAME [Factor]" or "NAME [Numeric]"'
    @params['backgroundExpression'] = 10
    @params['backgroundExpression', "description"] = "additive offset used in heatmaps"
    @params['transcriptTypes'] = ''
    @params['transcriptTypes', 'multi_selection'] = true
    @params['transcriptTypes', 'selected'] = 0

    @params['runGO'] = ['true', 'false']
    @params['runGO', 'description'] = "perform ORA and GSEA with Gene Ontology annotations"
    @params['pValThreshGO'] = 0.01
    @params['pValThreshGO', 'description'] = "pValue cut-off for ORA candidate gene selection"
    @params['log2RatioThreshGO'] = 0
    @params['log2RatioThreshGO', 'description'] = "log2 FoldChange cut-off for ORA candidate gene selection"
    @params['fdrThreshORA'] = 0.05
    @params['fdrThreshORA', 'description'] = "adjusted pValue cut-off for GO terms in ORA"
    @params['fdrThreshGSEA'] = 0.05
    @params['fdrThreshGSEA', 'description'] = "adjusted pValue cut-off for GO terms in GSEA"

    @params['specialOptions'] = ''
    @params['expressionName'] = ''
    @params['mail'] = ""
    @params['Rversion'] = ["Dev/R/4.3.0", "Dev/R/4.2.2", "Dev/R/4.2.0", "Dev/R/4.1.2", "Dev/R/4.1.0", "Dev/R/4.0.4", "Dev/R/4.0.3", "Dev/R/4.0.1"]
    @modules = ["Tools/samtools"]
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
     'Live Report [Link]'=>"http://fgcz-shiny.uzh.ch/exploreDEG?data=#{report_file}",
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
    command << run_RApp("EzAppDeseq2")
  end
end

if __FILE__ == $0

end
