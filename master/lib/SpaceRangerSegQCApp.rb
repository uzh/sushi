#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class SpaceRangerSegQCApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'SpaceRangerSegmentation'
    @analysis_category = 'Spatial'
    @description =<<-EOS
This wrapper runs <a href='https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/segmented-outputs',>space ranger segmentation</a> in Single-library analysis mode.
    EOS
    @required_columns = ['Name','Image']
    @required_params = ['name']
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '60'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '200'
    @params['scratch', "context"] = "slurm"
    @params['name'] = 'SpaceRangerSegmentation'
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'specify the commandline options for SpaceRanger; do not specify any option that is already covered by the dedicated input fields'
    @params['cmdOptions', "context"] = "SpaceRangerSegQC"
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @params['SpaceRangerVersion'] = ["Aligner/SpaceRanger/4.0.1"]
    @modules = ["Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def set_default_parameters
  end
  def next_dataset
    report_dir = File.join(@result_dir,"#{@dataset['Name']}")
    dataset = {
        'Name'=>@dataset['Name'],
        'ResultDir [File]'=>report_dir,
        'Report [Link]'=>File.join(report_dir, 'web_summary.html')       
      }.merge(extract_columns(@inherit_tags))
    dataset
  end
  def commands
    command = "module load  #{@params["SpaceRangerVersion"]}\n"
    command << run_RApp("EzAppSpaceRangerSegQC")
  end
end

if __FILE__ == $0

end
