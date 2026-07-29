#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class XeniumSeuratApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'XeniumSeurat'
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
    @params['name'] = 'XeniumSeurat'
    @params['minCounts'] = '10'
    @params['minCounts', 'description'] = 'Minimum counts per cell for QC filtering'
    @params['minFeatures'] = '5'
    @params['minFeatures', 'description'] = 'Minimum features per cell for QC filtering'
    @params['qcNmads'] = '3'
    @params['qcNmads', 'description'] = 'Number of MADs for QC outlier flagging (cell area, nucleus:cell ratio, etc.). Outliers are flagged and reported, not removed.'
    @params['clusterResolution'] = '0.5'
    @params['clusterResolution', 'description'] = 'Resolution for Seurat clustering (higher = more clusters)'
    @params['lambda'] = '0.8'
    @params['lambda', 'description'] = 'BANKSY spatial weighting (0-1). Higher values give more weight to spatial neighbors.'
    @params['nicheResolution'] = '0.5'
    @params['nicheResolution', 'description'] = 'Resolution for BANKSY spatial niche clustering (higher = more niches)'
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
      # Uveal Melanoma project references (human + melanoma cells) - pre-built RCTD format
      'p36005_UM_references/brain/brain_rctd.rds (human brain + melanoma)',
      'p36005_UM_references/breast/breast_rctd.rds (human breast + melanoma)',
      'p36005_UM_references/duodenum/duodenum_rctd.rds (human duodenum + melanoma)',
      'p36005_UM_references/eye/eye_rctd.rds (human eye + melanoma)',
      'p36005_UM_references/liver/liver_rctd.rds (human liver + melanoma)',
      'p36005_UM_references/skin/skin_rctd.rds (human skin + melanoma)',
      'p36005_UM_references/thyroid/thyroid_rctd.rds (human thyroid + melanoma)',
      # CTCL project references (human skin + CTCL malignant T cells TCM1/TCM2) - pre-built RCTD format
      'p28409_CTCL_references/skin_Ganier_TCM/skin_rctd.rds (human skin CTCL with malignant TCM1 TCM2)',
      'p28409_CTCL_references/skin_Ganier_noTCM/skin_rctd.rds (human skin no malignant T — for CBCL or non-CTCL samples)'
    ]
    @params['rctdReference', 'description'] = 'RCTD Reference atlas. Format: folder/file.rds (species tissue). WARNING: RCTD requires 200GB+ RAM.'
    @params['rctdFile'] = ''
    @params['rctdFile', 'description'] = 'Manual override: Full path to custom RCTD reference .rds file (leave empty to use dropdown selection)'
    @params['rctdUMImin'] = '10'
    @params['rctdUMImin', 'description'] = 'Minimum UMI count for RCTD annotation. Cells below this threshold will not be classified. Set to match minCounts (10) so low-RNA cells (T cells, neutrophils) that spillover hits hardest are retained.'
    @params['rctdClassFile'] = ''
    @params['rctdClassFile', 'description'] = 'Optional TSV (columns: cell_type, class) mapping reference cell types to coarse lineages. Sets RCTD class_df: maximises cell retention (same-class doublets become confident singlets) and unlocks splitMode=shift. Leave empty to disable.'
    @params['doSPLIT'] = false
    @params['doSPLIT', 'description'] = 'Run SPLIT spillover correction after RCTD (Bilous et al. 2026). Purifies transcript spillover/diffusion across cell boundaries; produces a before/after annotation panel. Requires an RCTD reference.'
    @params['splitMode'] = ['neighborhood', 'full', 'shift']
    @params['splitMode', 'description'] = 'SPLIT cell-selection. neighborhood (default, conservative): purify only if the secondary type is in the spatial neighbourhood. full: purify every contaminated cell (most aggressive). shift (experimental): transcriptomic-neighbourhood label correction, REQUIRES rctdClassFile.'
    @params['splitNeighborThreshold'] = '0.05'
    @params['splitNeighborThreshold', 'description'] = 'SPLIT neighborhood mode: minimum secondary-type weight in the spatial neighbourhood to trigger purification (balance_score_based threshold).'
    @params['coocRadius'] = '30'
    @params['coocRadius', 'description'] = 'Cell-type co-occurrence: spatial neighbour radius in microns for the neighbourhood-enrichment graph.'
    @params['coocNperm'] = '1000'
    @params['coocNperm', 'description'] = 'Cell-type co-occurrence: number of label permutations for the enrichment null (log2 enrichment + permutation p-value).'
    @params['rctdUMIminSigma'] = '300'
    @params['rctdUMIminSigma', 'description'] = 'Minimum UMI for the cells used to fit RCTD sigma (overdispersion), which drives every singlet/doublet/reject call. spacexr default 300, but Xenium medians run 100-300 so sigma gets fit on the count-rich tail. Lower it if log.txt warns that too few cells clear it.'
    @params['banksyDims'] = '12'
    @params['banksyDims', 'description'] = 'BANKSY: number of pca.banksy PCs used for niche clustering. Default 12 preserves historical behaviour; 30 (all PCs computed) uses the full embedding but shifts niche boundaries - change deliberately and note it when comparing against earlier runs.'
    @params['banksyKgeom'] = '30'
    @params['banksyKgeom', 'description'] = 'BANKSY: spatial neighbourhood size k_geom (the paper endorses 15-30).'
    @params['coocFdr'] = true
    @params['coocFdr', 'description'] = 'Cell-type co-occurrence: BH-correct the permutation p-values across cell-type pairs before starring the heatmap.'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    # Pin the R version, as every sibling spatial app does. An unversioned
    # "Dev/R" follows the Lmod default: it moved to 4.6.0 (Seurat 5.5.1), where
    # an upstream FindClusters change silently overwrote seurat_clusters, and
    # Dev/R/4.5.0 has no SPLIT package at all.
    # 4.6.0 only, deliberately. Dev/R/4.5.0 has no SPLIT package and carries a
    # stale XeniumSeurat.Rmd, so offering it would let a user pick a silently
    # degraded run. Add it back only once ezRun is deployed into that lib too.
    @params['Rversion'] = ["Dev/R/4.6.0"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def preprocess
    # Validate at SUBMIT time. Everything below used to fail (or silently
    # degrade) only inside the job, after SLURM had already granted 200 GB and
    # the pipeline had run for hours.
    if @params['rctdFile'].to_s != '' && !File.exist?(@params['rctdFile'].to_s)
      raise "rctdFile not found: #{@params['rctdFile']}"
    end
    if @params['rctdClassFile'].to_s != '' && !File.exist?(@params['rctdClassFile'].to_s)
      raise "rctdClassFile not found: #{@params['rctdClassFile']}"
    end
    if @params['doSPLIT']
      if @params['rctdFile'].to_s == '' &&
         (@params['rctdReference'].to_s == '' || @params['rctdReference'].to_s == 'None')
        raise "doSPLIT requires an RCTD reference: set rctdReference or rctdFile."
      end
      # splitMode 'shift' needs the class hierarchy. Without rctdClassFile,
      # RCTD emits no first_type_class/second_type_class and SPLIT's label-swap
      # dies inside dplyr - a failure the app's tryCatch used to swallow.
      if @params['splitMode'].to_s == 'shift' && @params['rctdClassFile'].to_s == ''
        raise "splitMode='shift' requires rctdClassFile (TSV: cell_type, class). Use splitMode='neighborhood' or supply the file."
      end
    end
  end
  def next_dataset
    # In SAMPLE mode, @dataset is a Hash (not Array) containing the current sample
    # During check phase, @dataset may be empty - return placeholder structure
    sample_name = @dataset['Name'] || 'placeholder'
    report_dir = File.join(@result_dir, sample_name)
    {'Name'=>sample_name,
     'Species'=>@dataset['Species'],
     'XeniumSeurat [File]'=>report_dir,
     'Static Report [Link]'=>File.join(report_dir, '00index.html'),
     # Expose the object itself so downstream apps and exploreSC links do not
     # have to guess the path (ScSeurat emits the equivalent 'SC Seurat [Link]').
     'SC Seurat [Link]'=>File.join(report_dir, 'scData.qs2')
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    command = "module load #{@params["Rversion"]}\n"
    command << run_RApp("EzAppXeniumSeurat")
  end
end
