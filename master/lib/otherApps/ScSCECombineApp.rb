#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class ScSCECombineApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'ScSCECombine'
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
    @params['node'] = ''
    @params['process_mode'] = 'DATASET'
    @params['name'] = 'SCReportMultipleSamplesSCE'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['tissue'] = []
    @params['tissue','multi_selection'] = true
    @params['tissue','all_selected'] = true
    @params['tissue', 'multi_selection_size'] = 10
    tissue = {}
    CSV.foreach("/srv/GT/databases/scGeneSets/all_cell_markers.txt", headers: true, col_sep: "\t") do |e|
      tissue[e["tissueType"]] = true
    end
    @params['tissue'] = tissue.keys.sort
    @params['tissue', 'description'] = 'Tissue the cells come from. Used in cell types identification for Human and Mouse organisms.'
    @params['resolution'] = '20'
    @params['resolution', 'description'] = 'A higher value will lead to more dense and less clusters.'
    @params['batchCorrection'] = 'true'
    @params['batchCorrection', 'description'] = "Perform batch correction?"
    @params['vars.regress'] = ['none', 'CellCycle']
    @params['vars.regress', 'description'] = 'Variable to regress when modelling the gene variance. Choose CellCycle if it is a bias.'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @params['Rversion'] = ["Dev/R/4.1.2", "Dev/R/4.1.0", "Dev/R/4.0.4", "Dev/R/4.0.3"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@params['name'],
     'Species'=>(dataset = @dataset.first and dataset['Species']),
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
    command = "module load #{@params["Rversion"]}\n"
    command << run_RApp("EzAppScSCECombine")
  end
end

if __FILE__ == $0

end
