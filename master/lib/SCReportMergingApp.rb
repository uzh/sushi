#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class SCReportMergingApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'SCReportMerging'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
    The report of merged single cell samples/plates<br/>
    EOS
    @required_columns = ['Name', 'Species', 'refBuild', 'refFeatureFile', 'Static Report']
    @required_params = []
    # optional params
    @params['cores'] = '4'
    @params['ram'] = '30'
    @params['scratch'] = '50'
    @params['name'] = 'SCReportMerging'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['scProtocol'] = ['10x', 'smart-Seq2']
    @params['pcs'] = 20
    @params['pcs', 'description'] = 'The maximal dimensions to use for reduction.'
    @params['x.low.cutoff'] = 0.0125
    @params['x.low.cutoff', 'description'] = 'Bottom cutoff on x-axis for identifying variable genes. 0.1 for smart-Seq2.'
    @params['x.high.cutoff'] = 3
    @params['x.high.cutoff', 'description'] = 'Top cutoff on x-axis for identifying variable genes. 8 for smart-Seq2.'
    @params['y.cutoff'] = 0.5
    @params['y.cutoff', 'description'] = 'Bottom cutoff on y-axis for identifying variable genes. 1 for msart-Seq2.'
    @params['resolution'] = 0.6
    @params['resolution', 'description'] = 'Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.'
    @params['batchCorrection'] = ['None', 'CCA']
    @params['batchCorrection', 'description'] = 'The batch correctio method: CCA or no correction'
    @params['cc'] = 20
    @params['cc', 'description'] = 'The number of canonical correlated subspaces.'
    @params['resolutionCCA'] = 0.6
    @params['resolutionCCA', 'description'] = 'Value of the resolution parameter for CCA correction version, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.'
    @params['chosenClusters'] =''
    @params['chosenClusters', 'description'] = 'The clusters to choose from each sample.In the format of sample1=cluster1,cluster2;sample2=cluster1,cluster2.'
    @params['all2allMarkers'] = false
    @params['all2allMarkers', 'description'] = 'Run all against all cluster comparisons?'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Dev/R", "Dev/Python/3.6.8"]
  end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@params['name'],
     'Species'=>(dataset = @dataset.first and dataset['Species']),
     'refBuild'=>@params['refBuild'],
     'refFeatureFile'=>@params['refFeatureFile'],
     'Static Report [Link]'=>report_link,
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
    run_RApp("EzAppSCReportMerging")
  end
end

if __FILE__ == $0

end
