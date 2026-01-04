#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class SeuratXeniumApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'SeuratXenium'
    @params['process_mode'] = 'SAMPLE'
    @analysis_category = 'Spatial'
    @description =<<-EOS
Seurat Analysis for Xenium Spatial Transcriptomics.<br/>
Includes QC, Normalization, Clustering, and RCTD Annotation.
    EOS
    @required_columns = ['Name','XeniumPath']
    @required_params = ['name']
    @params['cores'] = '8'
    @params['ram'] = '200'
    @params['scratch'] = '200'
    @params['name'] = 'SeuratXenium'
    @params['minCounts'] = '10'
    @params['minCounts', 'description'] = 'Minimum counts per cell for QC filtering'
    @params['minFeatures'] = '5'
    @params['minFeatures', 'description'] = 'Minimum features per cell for QC filtering'
    @params['Cluster_resolution'] = '0.5'
    @params['Cluster_resolution', 'description'] = 'Resolution for Seurat clustering (higher = more clusters)'
    @params['lambda'] = '0.8'
    @params['lambda', 'description'] = 'BANKSY spatial weighting (0-1). Higher values give more weight to spatial neighbors.'
    @params['Niche_resolution'] = '0.5'
    @params['Niche_resolution', 'description'] = 'Resolution for BANKSY spatial niche clustering (higher = more niches)'
    @params['rctdReference'] = ['None',
      # Mouse references
      'allen/allen_cortex_rctd.rds (mouse brain)',
      'azimuth/pan_mouse_pansci_rctd.rds (mouse pan-tissue)',
      'tabula_muris_senis/Bladder_rctd.rds (mouse)',
      'tabula_muris_senis/Bone_marrow_rctd.rds (mouse)',
      'tabula_muris_senis/Heart_rctd.rds (mouse)',
      'tabula_muris_senis/Kidney_rctd.rds (mouse)',
      'tabula_muris_senis/Large_intestine_rctd.rds (mouse)',
      'tabula_muris_senis/Limb_muscle_rctd.rds (mouse)',
      'tabula_muris_senis/Liver_rctd.rds (mouse)',
      'tabula_muris_senis/Lung_rctd.rds (mouse)',
      'tabula_muris_senis/Mammary_gland_rctd.rds (mouse)',
      'tabula_muris_senis/Pancreas_rctd.rds (mouse)',
      'tabula_muris_senis/Skin_rctd.rds (mouse)',
      'tabula_muris_senis/Spleen_rctd.rds (mouse)',
      'tabula_muris_senis/Thymus_rctd.rds (mouse)',
      'tabula_muris_senis/Tongue_rctd.rds (mouse)',
      'tabula_muris_senis/Trachea_rctd.rds (mouse)',
      # Human references - archmap
      'archmap/Glioblastoma_rctd.rds (human brain tumor)',
      'archmap/HLCA_rctd.rds (human lung)',
      'archmap/NSCLC_rctd.rds (human lung cancer)',
      # Human references - celltypist
      'celltypist/Blood_rctd.rds (human)',
      'celltypist/Bone_marrow_rctd.rds (human)',
      'celltypist/Heart_rctd.rds (human)',
      'celltypist/Hippocampus_rctd.rds (human brain)',
      'celltypist/Intestine_rctd.rds (human)',
      'celltypist/Kidney_rctd.rds (human)',
      'celltypist/Liver_rctd.rds (human)',
      'celltypist/Lung_rctd.rds (human)',
      'celltypist/Lymph_node_rctd.rds (human)',
      'celltypist/Pancreas_rctd.rds (human)',
      'celltypist/Skeletal_muscle_rctd.rds (human)',
      'celltypist/Spleen_rctd.rds (human)',
      # Human references - disco (complete)
      'disco/adipose_rctd.rds (human)',
      'disco/adrenal_gland_rctd.rds (human)',
      'disco/bladder_rctd.rds (human)',
      'disco/blood_rctd.rds (human)',
      'disco/bone_marrow_rctd.rds (human)',
      'disco/brain_rctd.rds (human)',
      'disco/breast_rctd.rds (human)',
      'disco/eye_rctd.rds (human)',
      'disco/fallopian_tube_rctd.rds (human)',
      'disco/gingiva_rctd.rds (human)',
      'disco/heart_rctd.rds (human)',
      'disco/intestine_rctd.rds (human)',
      'disco/kidney_rctd.rds (human)',
      'disco/liver_cell_rctd.rds (human)',
      'disco/lung_rctd.rds (human)',
      'disco/ovary_rctd.rds (human)',
      'disco/pancreas_cell_rctd.rds (human)',
      'disco/placenta_rctd.rds (human)',
      'disco/skeletal_muscle_rctd.rds (human)',
      'disco/skin_rctd.rds (human)',
      'disco/stomach_rctd.rds (human)',
      'disco/testis_rctd.rds (human)',
      'disco/thymus_rctd.rds (human)',
      'disco/tonsil_rctd.rds (human)',
      # Human disease-specific references - disco
      'disco/AD_frontal_cortex_parenchyma_rctd.rds (human Alzheimer)',
      'disco/COVID-19_blood_rctd.rds (human COVID-19)',
      'disco/Crohns_disease_ileum_rctd.rds (human Crohn)',
      'disco/PDAC_pancreas_rctd.rds (human pancreatic cancer)',
      'disco/type_1_diabetes_pancreas_rctd.rds (human T1D)',
      'disco/type_2_diabetes_pancreas_rctd.rds (human T2D)',
      # Uveal Melanoma project references (human + melanoma cells)
      'p36005_UM_references/brain/brain_ref_downsampled.qs2 (human brain + melanoma)',
      'p36005_UM_references/breast/breast_melanoma_reference_raw_counts.qs2 (human breast + melanoma)',
      'p36005_UM_references/duodenum/duod_ref_downsampled.qs2 (human duodenum + melanoma)',
      'p36005_UM_references/eye/eye_ref_downsampled.qs2 (human eye/retina + melanoma)',
      'p36005_UM_references/liver/ref_downsampled.qs (human liver + melanoma)',
      'p36005_UM_references/skin/skin_ref_downsampled.qs2 (human skin + melanoma)',
      'p36005_UM_references/thyroid/thyroid_ref_downsampled.qs2 (human thyroid + melanoma)'
    ]
    @params['rctdReference', 'description'] = 'RCTD Reference atlas. Format: folder/file.rds (species tissue). WARNING: RCTD requires 200GB+ RAM.'
    @params['rctdFile'] = ''
    @params['rctdFile', 'description'] = 'Manual override: Full path to custom RCTD reference .rds file (leave empty to use dropdown selection)'
    @params['rctdUMImin'] = '20'
    @params['rctdUMImin', 'description'] = 'Minimum UMI count for RCTD annotation. Cells below this threshold will not be classified.'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def next_dataset
    # In SAMPLE mode, @dataset is a Hash (not Array) containing the current sample
    # During check phase, @dataset may be empty - return placeholder structure
    sample_name = @dataset['Name'] || 'placeholder'
    report_dir = File.join(@result_dir, sample_name)
    {'Name'=>sample_name,
     'ReportData [File]'=>report_dir,
     'Report [Link]'=>File.join(report_dir, '00index.html'),
     'Species'=>@dataset['Species']
    }
  end
  def commands
    run_RApp("EzAppSeuratXenium")
  end
end
