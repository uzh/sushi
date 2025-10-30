#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class CellBenderApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'CellBender'
    @params['process_mode'] = 'SAMPLE'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
CellBender is a software package for eliminating technical artifacts from high-throughput single-cell RNA sequencing (scRNA-seq) data. It is often also referred to by its French name, Le plieur de cellules. ToolLink: https://github.com/broadinstitute/CellBender<br/>
    EOS
    @required_columns = ['Name', 'Species', 'refBuild', 'refFeatureFile', 'CountMatrix']
    @required_params = ['name']
    # optional params
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '30'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '100'
    @params['scratch', "context"] = "slurm"
    @params['gpu'] = '1'
    @params['gpu', "context"] = "CellBender"
    @params['name'] = 'CellBender'
    @params['gpuMode'] = false
    @params['gpuMode', "context"] = "CellBender"
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'specify other commandline options for CellBender'
    @params['cmdOptions', "context"] = "CellBender"
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "referfence genome assembly"
    @params['refFeatureFile'] = 'genes.gtf'
    @params['refFeatureFile', "context"] = "CellBender"
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def next_dataset
    report_file = File.join(@result_dir, "#{@dataset['Name']}")
    report_link = File.join(report_file, 'cellbender_report.html')
    {'Name'=>@dataset['Name'],
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'refFeatureFile'=>@params['refFeatureFile'],
     'Static Report [Link]'=>report_link,
     'ResultDir [File]'=>report_file,
     'CountMatrix [Link]'=>File.join(report_file,'cellbender_filtered_seurat.h5'),
     'UnfilteredCountMatrix [Link]'=>File.join(report_file,'cellbender_raw_seurat.h5'),
    }.merge(extract_columns(@inherit_tags))
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('refFeatureFile')
      @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    end
  end
  def commands
    run_RApp("EzAppCellBender",conda_env: "gi_cellbender_0.3.2")
  end
end

if __FILE__ == $0

end
