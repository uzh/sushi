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
    An ENACT app for VisiumHD data.<br/>
    The Pipeline is adapted from https://github.com/Sanofi-Public/enact-pipeline/tree/main<br/>
    ENACT is a tool for converting binned spatial transcriptomics data into single-cell format.<br/>
    Their paper can be found here: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btaf094/8063614<br/>
    This app is running fully on Python, therefore, check the implementation in ezPyz: https://github.com/fgcz/EzPyzApps . There you can find the required information and the implementation of the app.<br/>
    The compiled code lives in the conda environment `gi_enact`<br/>
    EOS
    @required_columns = ['Name','BinnedOutputs2um','SourceImage']
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

    @params['destripe_norm'] = true
    @params['destripe_norm', 'description'] = 'Run the destriping step (same as in bin2cell)'
    @params['prob_thresh'] = 0.005
    @params['prob_thresh', 'description'] = 'Probability threshold for stardist cell detection'
    @params['expand_by_nbins'] = 2
    @params['expand_by_nbins', 'description'] = 'Distance (in bin lengths) to expand the nuclei to obtain cell areas (0 for false)'
    @params['n_hvg'] = 1000
    @params['n_hvg', 'description'] = 'Number of highly variable genes to use (20k means all genes)'

    @params['bin_to_cell_method'] = ['weighted_by_area', 'naive', 'weighted_by_cluster', 'weighted_by_gene']
    @params['bin_to_cell_method', 'description'] = 'Method to assign bins to cells'
    @params['bin_to_cell_method', 'default'] = 'weighted_by_area'
    @params['bin_to_cell_method', 'single_selection'] = true
    @params['cell_annotation_method'] = ['celltypist', 'cellassign']
    @params['cell_annotation_method','description'] = 'Use the cell-wise transcript counts to infer the cell labels/ phenotypes using methods used for single-cell RNA seq analysis'    
    @params['cell_annotation_method', 'default'] = 'celltypist'
    @params['cell_annotation_method', 'single_selection'] = true
    @params['cell_typist_model'] = ['Immune_All_Low.pkl', 
                                    'Immune_All_High.pkl', 
                                    'Adult_COVID19_PBMC.pkl', 
                                    'Adult_CynomolgusMacaque_Hippocampus.pkl', 
                                    'Adult_Human_MTG.pkl', 
                                    'Adult_Human_PancreaticIslet.pkl', 
                                    'Adult_Human_PrefrontalCortex.pkl', 
                                    'Adult_Human_Skin.pkl', 
                                    'Adult_Human_Vascular.pkl', 
                                    'Adult_Mouse_Gut.pkl', 
                                    'Adult_Mouse_OlfactoryBulb.pkl', 
                                    'Adult_Pig_Hippocampus.pkl', 
                                    'Adult_RhesusMacaque_Hippocampus.pkl', 
                                    'Autopsy_COVID19_Lung.pkl', 
                                    'COVID19_HumanChallenge_Blood.pkl', 
                                    'COVID19_Immune_Landscape.pkl', 
                                    'Cells_Adult_Breast.pkl', 
                                    'Cells_Fetal_Lung.pkl', 
                                    'Cells_Human_Tonsil.pkl', 
                                    'Cells_Intestinal_Tract.pkl', 
                                    'Cells_Lung_Airway.pkl', 
                                    'Developing_Human_Brain.pkl', 
                                    'Developing_Human_Gonads.pkl', 
                                    'Developing_Human_Hippocampus.pkl', 
                                    'Developing_Human_Organs.pkl', 
                                    'Developing_Human_Thymus.pkl', 
                                    'Developing_Mouse_Brain.pkl', 
                                    'Developing_Mouse_Hippocampus.pkl', 
                                    'Fetal_Human_AdrenalGlands.pkl', 
                                    'Fetal_Human_Pancreas.pkl', 
                                    'Fetal_Human_Pituitary.pkl', 
                                    'Fetal_Human_Retina.pkl', 
                                    'Fetal_Human_Skin.pkl', 
                                    'Healthy_Adult_Heart.pkl', 
                                    'Healthy_COVID19_PBMC.pkl', 
                                    'Healthy_Human_Liver.pkl', 
                                    'Healthy_Mouse_Liver.pkl', 
                                    'Human_AdultAged_Hippocampus.pkl', 
                                    'Human_Colorectal_Cancer.pkl', 
                                    'Human_Developmental_Retina.pkl', 
                                    'Human_Embryonic_YolkSac.pkl', 
                                    'Human_Endometrium_Atlas.pkl', 
                                    'Human_IPF_Lung.pkl', 
                                    'Human_Longitudinal_Hippocampus.pkl', 
                                    'Human_Lung_Atlas.pkl', 
                                    'Human_PF_Lung.pkl', 
                                    'Human_Placenta_Decidua.pkl', 
                                    'Lethal_COVID19_Lung.pkl', 
                                    'Mouse_Dentate_Gyrus.pkl', 
                                    'Mouse_Isocortex_Hippocampus.pkl', 
                                    'Mouse_Postnatal_DentateGyrus.pkl', 
                                    'Mouse_Whole_Brain.pkl', 
                                    'Nuclei_Lung_Airway.pkl', 
                                    'Pan_Fetal_Human.pkl'
                                  ]
    @params['cell_typist_model', 'description'] = 'Specify the model to use with CellTypist (if selected above), models can be found here: https://www.celltypist.org/models'
    @modules = []
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def next_dataset    
    dir_name = "#{@params['name']}_#{@dataset['Name']}"
    report_dir = File.join(@result_dir, dir_name) #might need to add equiv of ENACT_output_files/ENACTApp_VisiumHD_Colon_Cancer
    {'Name'=> dir_name,
    'ENACT [File]'=> report_dir,
    'Anndata [Link]' => File.join(report_dir, "chunks", "#{@params['bin_to_cell_method']}", "#{@params['cell_annotation_method']}_results", 'cells_adata.h5'),
    'TissUUmap [Link]'=> File.join(report_dir, 'tmap')
    }
  end
  def commands
    run_PyApp("ENACT",conda_env: 'gi_enact')
  end
end
