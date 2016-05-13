#!/usr/bin/env ruby
# encoding: utf-8
Version = '20160513'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class DEXSeqApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'DEXSeqApp'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Alternative_Splicing'
    @description =<<-EOS
    Inference of differential exon usage in RNA-Seq<br/>
<a href='https://bioconductor.org/packages/release/bioc/html/DEXSeq.html'>DEXSeq</a><br/>
    EOS
    @required_columns = ['Name','BAM','BAI', 'Species','refBuild', 'refFeatureFile']
    @required_params = ['grouping', 'sampleGroup', 'refGroup','refBuild']
    # optional params
    @params['cores'] = '1'
    @params['ram'] = '2'
    @params['scratch'] = '10'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['grouping'] = '' 
    @params['sampleGroup'] = '' 
    @params['refGroup'] = '' 
    @params['expressionName'] = ''
    @params['specialOptions'] = ''
    @params['mail'] = ""
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
     'Report [File]'=>report_file,
     'Html [Link]'=>report_link,
    }
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('refFeatureFile')
      @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    end
  end
  def commands
    run_RApp("EzAppDEXSeqAnalysis")
  end
end

if __FILE__ == $0

end

