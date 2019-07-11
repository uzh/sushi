#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class SCMultipleSamplesAndGroupsApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'SCMultipleSamplesAndGroups'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
    The report of merged single cell samples/plates<br/>
    EOS
    @required_columns = ['Name', 'Species', 'refBuild', 'refFeatureFile', 'Static Report']
    @required_params = []
    # optional params
    @param['cores'] = '4'
    @param['ram'] = '50'
    @param['scratch'] = '50'
    @param['node'] = ''
    @param['process_mode'] = 'DATASET'
    @param['samples'] = 'Epilepsy_13122018,GBM_28_08_2018'
    @param['name'] = 'SCReportMerging'
    @param['refBuild'] = ref_selector
    @param['refFeatureFile'] = 'genes.gtf'
    @param['npcs'] = '20'
    @param['resolution'] = '0.6'
    @param['batchCorrection'] = 'true'
    @param['chosenClusters'] = ''
    @param['all2allMarkers'] = 'false'
    @param['specialOptions'] = ''
    @param['mail'] = ""
    @param['Rversion'] = 'Dev/R/3.6.0'
    @modules = ["Dev/R", "Dev/Python/3.6.8"]
  end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@@params['name'],
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
    run_RApp("EzAppSCMultipleSamplesAndGroups")
  end
end

if __FILE__ == $0

end
