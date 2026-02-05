#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class VelocytoApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'Velocyto'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
This wrapper runs <a href='http://velocyto.org/',>velocyto</a> in Single-library analysis mode.
    EOS
    @required_columns = ['Name', 'ResultDir', 'AlignmentFile']
    @required_params = ['name', 'refBuild']
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '30'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '200'
    @params['scratch', "context"] = "slurm"
    @params['name'] = 'Velocyto'
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "referfence genome assembly"
    @params['refFeatureFile'] = 'genes.gtf'
    @params['refFeatureFile', "context"] = "Velocyto"
    @params['featureLevel'] = 'genes'
    @params['featureLevel', "context"] = "Velocyto"
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'specify the commandline options for Velocyto  do not specify any option that is already covered by the dedicated input fields'
    @params['cmdOptions', "context"] = "Velocyto"
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Dev/R", "Tools/samtools"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    @params['featureLevel'] = @dataset[0]['featureLevel']
  end
  def next_dataset
      dataset = {
        'Name'=>@dataset['Name'],
        'LOOM [File]'=>File.join(@result_dir, "#{@dataset['Name']}.loom"),
        'Species'=>@dataset['Species'],
        'refBuild'=>@params['refBuild'],
        'refFeatureFile'=>@params['refFeatureFile'],
        'featureLevel'=>@params['featureLevel'],
        'Read Count'=>@dataset['Read Count']
      }.merge(extract_columns(@inherit_tags))
    dataset
  end
  def commands
    command = "module load  #{@params["CellRangerVersion"]}\n"
    command << run_RApp("EzAppVeloCyto", conda_env: "gi_velocyto")
  end
end

if __FILE__ == $0

end
