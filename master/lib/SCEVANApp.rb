#!/usr/bin/env ruby
# encoding: utf-8

################################################################################
# FGCZ SUSHI App: SCEVAN CNV Analysis
# Ruby on Rails frontend definition for SCEVAN pipeline
#
# IMPORTANT: RESOURCE REQUIREMENTS
#
# SCEVAN is a highly memory-intensive tool requiring substantial computational resources:
#   - RAM: Minimum 500 GB, recommended 600-800 GB for typical scRNA datasets
#   - Cores: Minimum 48, can utilize up to 64 cores efficiently
#   - Runtime: Can take 4-12 hours for datasets with 50K-200K cells
#   - Scratch: 200 GB for intermediate files and SCEVAN outputs
#
# Ensure adequate resources are allocated before submission!
# Consider testing with a small subset first.
################################################################################

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class SCEVANApp <  SushiFabric::SushiApp
  def initialize
    super

    # ========================================================================
    # BASIC APP CONFIGURATION
    # ========================================================================

    # App name (must match EzAppSCEVANApp in R)
    @name = 'SCEVANApp'

    # Process mode: DATASET = process all samples together
    @params['process_mode'] = 'DATASET'

    # Analysis category
    @analysis_category = 'SingleCell'

    # Description shown in SUSHI web interface
    @description =<<-EOS
<b>SCEVAN (Single CEll Variational Aneuploidy aNalysis)</b> infers copy number variations (CNV) from single-cell RNA-seq data.<br/>
<br/>
<b>Key features:</b><br/>
<ul>
  <li>Classifies cells as malignant vs. non-malignant based on CNV patterns</li>
  <li>Identifies tumor subclones with distinct copy number architectures</li>
  <li>Generates chromosome-level CNV profiles for each cell</li>
  <li>Supports both human and mouse organisms</li>
  <li>Requires user-specified reference (normal) cell types for comparison</li>
</ul>
<br/>
<b>⚠️ Resource Requirements:</b><br/>
<ul>
  <li><b>RAM:</b> 500-800 GB (SCEVAN is memory-intensive!)</li>
  <li><b>Cores:</b> 48-64 recommended</li>
  <li><b>Runtime:</b> 4-12 hours for typical datasets</li>
</ul>
<br/>
<b>Input:</b> Seurat object (qs2 or RDS format) with cell type annotations<br/>
<b>Output:</b> HTML report with CNV classifications, UMAPs, heatmaps, and results tables<br/>
<br/>
For more information, see: <a href='https://github.com/AntonioDeFalco/SCEVAN'>SCEVAN GitHub</a>
    EOS

    # ========================================================================
    # REQUIRED DATASET COLUMNS
    # ========================================================================
    @required_columns = ['Name', 'SeuratObject', 'Species']

    # ========================================================================
    # REQUIRED PARAMETERS
    # ========================================================================
    @required_params = ['name', 'cellTypeColumn', 'referenceCellTypes']

    # ========================================================================
    # COMPUTATIONAL RESOURCES
    # ========================================================================
    # Number of CPU cores (SCEVAN benefits from high parallelization)
    @params['cores'] = '48'

    # RAM in GB (SCEVAN requires >500 GB for typical datasets)
    @params['ram'] = ['500', '600', '800']

    # Scratch space in GB
    @params['scratch'] = '200'

    # ========================================================================
    # OUTPUT NAME
    # ========================================================================
    @params['name'] = 'SCEVAN_CNV_Analysis'

    # ========================================================================
    # EMAIL NOTIFICATION
    # ========================================================================
    @params['mail'] = ""

    # ========================================================================
    # SCEVAN ANALYSIS PARAMETERS
    # ========================================================================

    # --- CELL TYPE COLUMN ---
    # Column in Seurat metadata containing cell type annotations
    @params['cellTypeColumn'] = 'predicted.id'
    @params['cellTypeColumn', 'description'] = 'Metadata column containing cell type annotations. Common options: predicted.id, cellTypeIntegrated, Azimuth.celltype.l2, manual_cell_type_L1, seurat_clusters'

    # --- REFERENCE CELL TYPES ---
    # Cell types to use as normal reference (comma or newline separated)
    @params['referenceCellTypes'] = 'CD4-positive, alpha-beta T cell
CD8-positive, alpha-beta T cell
B cell
natural killer cell
plasma cell'
    @params['referenceCellTypes', 'description'] = 'Cell types to use as normal reference for CNV calling (comma or newline separated). Use exact names from your cell type column. Examples for human: CD4+ T cells, CD8+ T cells, B cells, NK cells, plasma cells'
    @params['referenceCellTypes', 'rows'] = 6

    # --- SAMPLE COLUMN ---
    # Column containing sample identifiers
    @params['sampleColumn'] = 'Sample'
    @params['sampleColumn', 'description'] = 'Metadata column containing sample identifiers (for multi-sample analysis)'

    # --- ORGANISM ---
    @params['organism'] = ['human', 'mouse']
    @params['organism', 'description'] = 'Organism: human or mouse (determines reference genome for CNV calling)'

    # --- ENABLE SUBCLONES ANALYSIS ---
    @params['subclones'] = true
    @params['subclones', 'description'] = 'Enable tumor subclone identification and phylogenetic analysis'

    # --- BETA VEGA PARAMETER ---
    @params['betaVega'] = '0.5'
    @params['betaVega', 'description'] = 'Segmentation granularity parameter (0-1). Lower values = more fine-grained segments. Default: 0.5'

    # ========================================================================
    # ENVIRONMENT MODULES
    # ========================================================================
    @modules = ["Dev/R"]

    # ========================================================================
    # CONDA ENVIRONMENT (optional)
    # ========================================================================
    # Uncomment if SCEVAN requires specific conda environment
    # @conda_env = 'scevan_env'

    # ========================================================================
    # INHERITED TAGS AND COLUMNS
    # ========================================================================
    @inherit_tags = ["Factor", "B-Fabric"]
  end

  # ==========================================================================
  # DEFINE OUTPUT DATASET STRUCTURE
  # ==========================================================================
  def next_dataset
    # Output directory path
    report_dir = File.join(@result_dir, @params['name'])

    # Define output dataset structure
    dataset = {
      'Name' => @params['name'],
      'Report [Link]' => File.join(report_dir, '00index.html'),
      'ReportData [File]' => report_dir,
      'Species' => (@dataset.first and @dataset.first['Species']),
      'SeuratObject [File]' => File.join(report_dir, 'seurat_object.qs2'),
      'SCEVANResults [File]' => File.join(report_dir, 'scevan_results_list.qs2')
    }

    return dataset
  end

  # ==========================================================================
  # GENERATE R COMMAND
  # ==========================================================================
  def commands
    run_RApp("EzAppSCEVANApp")
  end
end

################################################################################
# USAGE NOTES
################################################################################
#
# 1. CELL TYPE COLUMN EXAMPLES:
#    After ScSeuratCombinedLabelClusters:
#      - cellTypeIntegrated (primary output)
#      - ident
#
#    From Azimuth:
#      - predicted.id (most common)
#      - predicted.celltype.l1 (broad)
#      - predicted.celltype.l2 (fine-grained)
#
#    From SingleR:
#      - annot_bmm_broad (BoneMarrowMap)
#      - annot_immgen (ImmGen reference)
#
#    Manual/Other:
#      - manual_cell_type_L1, manual_cell_type_L2
#      - seurat_clusters (fallback)
#
# 2. REFERENCE CELL TYPES:
#    Must match EXACTLY the names in your cell type column
#    Common human examples:
#      - T cells: CD4-positive, alpha-beta T cell
#                 CD8-positive, alpha-beta T cell
#      - B cells: B cell, naive B cell, memory B cell, plasma cell
#      - Myeloid: monocyte, macrophage, dendritic cell
#      - Other: natural killer cell, erythrocyte
#
#    Mouse examples will depend on the reference used
#
# 3. TESTING:
#    Test with small datasets first (e.g., subset to 5000-10000 cells)
#    Monitor memory usage - SCEVAN can require >500 GB RAM
#
# 4. TROUBLESHOOTING:
#    - If job fails with memory error: increase RAM allocation
#    - If no normal cells found: check cell type names match exactly
#    - If <5% normal cells: consider adding more reference cell types
#
################################################################################
