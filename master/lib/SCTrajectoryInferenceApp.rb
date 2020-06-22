#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class SCTrajectoryInferenceApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'SCTrajectoryInference'
    @params['process_mode'] = 'SAMPLE'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
Trajectory inference analysis for single cell data<br/>
    EOS
    @required_columns = ['Name', 'refBuild', 'Report', 'ResultDir']
    @required_params = ['name']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '250'
    @params['name'] = 'SCTrajectoryInference'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['start_id'] = '0'
    @params['end_id'] = ''
    @params['start_n'] = '1'
    @params['end_n'] = '1'
    @params['TI_method'] = ''
    @params['diff_Branch'] = ''
    @params['diff_Branch_Point'] = ''
    @params['specialOptions'] = ''
    @params['mail'] = ''
    @modules = ["Dev/R", "Dev/Python/3.6.8"]
  end
  def next_dataset
    report_file = File.join(@result_dir, "#{@dataset['Name']}_SCTrajectoryInference")
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@dataset['Name'],
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
    run_RApp("EzAppTrajectoryInference")
  end
end

if __FILE__ == $0

end

