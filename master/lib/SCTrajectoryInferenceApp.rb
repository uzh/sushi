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
    @required_columns = ['Name', 'Report', 'refBuild', 'refFeatureFile']
    @required_params = ['name']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '250'
    @params['name'] = 'SCTrajectoryInference'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['start_id'] = '0'
    @params['start_id', 'description'] = 'Start cluster(s)'
    @params['end_id'] = 'none'
    @params['end_id', 'description'] = 'End cluster(s)'
    @params['start_n'] = '1'
    @params['start_n', 'description'] = 'The number of start states'
    @params['end_n'] = '1'
    @params['end_n', 'description'] = 'The number of end states'
    @params['TI_method'] = 'none'
    @params['TI_method', 'description'] = 'Trajectory inference method(s)'
    @params['show_genes'] = 'none'
    @params['show_genes', 'description'] = 'Genes to show along the trajectory. (For example: SCGB1B1,SCGB3A1)'
    @params['root_expression'] = 'none'
    @params['root_expression', 'description'] = 'Rooting the trajectory by using genes that are highly expressed at the start of it. (For example: SCGB1B1,SCGB3A1)'
    @params['diff_Branch'] = 'none'
    @params['diff_Branch', 'description'] = 'Method and branch name to extract dysregulated genes from. (For example: Slingshot,3)'
    @params['diff_Branch_Point'] = 'none'
    @params['diff_Branch_Point', 'description'] = 'Method and branching point name to extract dysregulated genes from. (For example: Slingshot,3)'
    @params['specialOptions'] = ''
    @params['mail'] = ''
    @modules = ["Dev/R", "Dev/Python"]
  end
  def next_dataset
    report_file = File.join(@result_dir, "#{@dataset['Name']}_SCTrajectoryInference")
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@dataset['Name'],
     'Report [File]'=>report_file,
     'Static Report [Link]'=>report_link
    }
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('refFeatureFile')
      @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    end
  end
  def commands
    run_RApp("EzAppSCTrajectoryInference")
  end
end

if __FILE__ == $0

end

