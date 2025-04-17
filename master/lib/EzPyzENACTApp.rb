#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class EzPyzENACTApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'EzPyzENACTApp'
    @params['process_mode'] = 'SAMPLE'
    @analysis_category = 'Spatial'
    @description =<<-EOS
    An ENACT app for VisiumHD data<br/>
    EOS
    @required_columns = ['Name','BinnedOutputs2um','SourceImage','SpaceRanger']
    @required_params = ['name']
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '100'
    @params['name'] = 'ENACTApp'
    @params['mail'] = ''    
    @params['segmentation'] = true
    @params['segmentation', 'description'] = 'Run the segmentation step'
    @params['bin_to_geodataframes'] = true
    @params['bin_to_geodataframes', 'description'] = 'Run the bin_to_geodataframes step'
    @params['bin_to_cell_assignment'] = true
    @params['bin_to_cell_assignment', 'description'] = 'Run the bin_to_cell_assignment step'
    @params['cell_type_annotation'] = true
    @params['cell_type_annotation', 'description'] = 'Run the cell_type_annotation step'
    @params['bin_to_cell_method'] = ['weighted_by_area', 'naive', 'weighted_by_cluster', 'weighted_by_gene']
    @params['bin_to_cell_method', 'description'] = 'Method to assign bins to cells'
    @params['bin_to_cell_method', 'default'] = 'weighted_by_area'
    @params['bin_to_cell_method', 'single_selection'] = true
    @params['cell_annotation_method'] = ['celltypist', 'cellassign']
    @params['cell_annotation_method','description'] = 'Use the cell-wise transcript counts to infer the cell labels/ phenotypes using methods used for single-cell RNA seq analysis'    
    @params['cell_annotation_method', 'default'] = 'celltypist'
    @params['cell_annotation_method', 'single_selection'] = true
    @params['cell_typist_model'] = ''
    @params['cell_typist_model', 'description'] = 'Specify the model to use with CellTypist (if selected above)'
    @modules = []
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def next_dataset    
    dir_name = "#{@params['name']}_#{@dataset['Name']}"
    report_dir = File.join(@result_dir, dir_name) #might need to add equiv of ENACT_output_files/ENACTApp_VisiumHD_Colon_Cancer
    {'Name'=> dir_name,
    'ENACT [File]'=> report_dir,
    'Anndata [File]' => File.join(report_dir, "#{@params['bin_to_cell_method']}", "#{@params['cell_annotation_method']}_results", 'cells_adata.csv'),
    'TissUUmap [Link]'=> File.join(report_dir, 'tmap')
    }
  end
  def commands
    run_PyApp("ENACT",conda_env: 'tmp_enact')
  end
end
