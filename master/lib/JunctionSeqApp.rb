#!/usr/bin/env ruby
# encoding: utf-8
Version = '20180412-100423'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class JunctionSeqApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'JunctionSeqApp'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Alternative_Splicing'
    @description =<<-EOS
    Inference of differential exon or Splice-Junction usage in RNA-Seq<br/>
<a href='https://bioconductor.org/packages/release/bioc/html/JunctionSeq.html'>JunctionSeq</a><br/>
    EOS
    @required_columns = ['Name','BAM','BAI', 'Species','refBuild', 'refFeatureFile']
    @required_params = ['grouping', 'sampleGroup', 'refGroup','refBuild','paired', 'strandMode']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '50'
    @params['scratch'] = '100'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['paired'] = false
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['grouping'] = ''
    @params['sampleGroup'] = ''
    @params['sampleGroup', 'description'] = 'sampleGroup should be different from refGroup'
    @params['refGroup'] = ''
    @params['refGroup', 'description'] = 'refGroup should be different from sampleGroup'
    @params['fdr'] = '0.05'
    @params['maxReadLength'] = '151'
    @params['minCount'] = '6'
    @params['maxGenes'] = '1000'
    @params['refGroup', 'description'] = 'maximum number of genes for plotting'
    @params['expressionName'] = ''
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Dev/jdk", "Dev/R/4.1.2"]
  end
  def next_dataset
    @comparison = "#{@params['sampleGroup']}--over--#{@params['refGroup']}"
    @params['comparison'] = @comparison
    @params['name'] = @comparison
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, 'finalOutput/testForDU.html')
    {'Name'=>@comparison,
     'Species'=>(dataset = @dataset.first and dataset['Species']),
     'refBuild'=>@params['refBuild'],
     'Report [File]'=>report_file,
     'Html [Link]'=>report_link,
    }
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('refFeatureFile')
      @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    end
    if dataset_has_column?('paired')
      @params['paired'] = @dataset[0]['paired']
    end
    if dataset_has_column?('strandMode')
      @params['strandMode'] = @dataset[0]['strandMode']
    end
  end
  def commands
    run_RApp("EzAppJunctionSeq")
  end
end

if __FILE__ == $0

end
