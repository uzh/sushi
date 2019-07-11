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
    @params['cores'] = '4'
    @params['ram'] = '50'
    @params['scratch'] = '50'
    @params['node'] = ''
    @params['process_mode'] = 'DATASET'
    @params['name'] = 'SCReportMerging'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['npcs'] = '20'
    @params['resolution'] = '0.6'
    @params['batchCorrection'] = 'true'
    @params['chosenClusters'] = ''
    @params['all2allMarkers'] = 'false'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @params['Rversion'] = 'Dev/R/3.6.0'
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
    run_RApp("EzAppSCMultipleSamplesAndGroups")
  end
end

if __FILE__ == $0

end
