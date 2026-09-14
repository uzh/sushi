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
    @required_columns = ['Name','Count', 'Species', 'libName']
    @required_params = []
    # optional params
    @params['cores'] = ['1', '2']
    @params['cores', "context"] = "slurm"
    @params['ram'] = ['4', '7', '50']
    @params['ram', "context"] = "slurm"
    @params['scratch'] = ['10', '20']
    @params['scratch', "context"] = "slurm"
    @params['name'] = 'ExploreMageckCounts'
    @params['normMethod'] = ['deseq2', 'tmm', 'cpm', 'logMean']
    @params['normMethod', "description"] = "count normalisation: deseq2 (size factors, as MAGeCK), edgeR tmm/cpm, or logMean"
    @params['refGroup'] = ''
    @params['refGroup', "description"] = "reference/plasmid/T0 condition used as baseline for essential-gene depletion + ROC; leave empty to skip that analysis"
    @params['backgroundExpression'] = 5
    @params['backgroundExpression', "description"] = "pseudo-count added before the log2 transform"
    @params['topGeneSize'] = 100
    @params['topGeneSize', "description"] = "number of most-variable sgRNAs shown in the top-feature heatmap"
    @params['nSampleClusters'] = 6
    @params['nSampleClusters', "description"] = "number of sample clusters annotated in the high-variance heatmap"
    @params['mail'] = ""
    @modules = ["Dev/R"]
    @inherit_columns = ["Order Id"]
  end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@params['name'],
     'Species'=>(dataset = @dataset.first and dataset['Species']),
     'libName'=>(dataset = @dataset.first and dataset['libName']),
     'Static Report [Link]'=>report_link,
     'Report [File]'=>report_file,
    }.merge(extract_columns(colnames: @inherit_columns))
  end
  def commands
    run_RApp("EzAppExploreMageckCounts")
  end
end

if __FILE__ == $0

end

