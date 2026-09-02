#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class ScSeuratCombineApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'ScSeuratCombine'
    # Seurat (v5) unconditional. CCA/RPCA vs Harmony integration mutually
    # exclusive per param$integrationMethod. Enrichr gated on enrichrDatabase;
    # decoupleR/dorothea/progeny gated on computePathwayTFActivity.
    # NOTE: unlike ScSeuratApp, this app does NOT declare params for per-sample QC
    # (doublets/ambient RNA/emptyDrops) or most cell-type annotation tools
    # (AUCell/SingleR/sc-type/mLLMCelltype/CyteTypeR/Azimuth) -- those ran upstream
    # in ScSeuratApp already. It only reads a cached cellxgeneResults.rds if
    # present, never computes it, so STACAS/schard aren't invoked here either.
    @citation = [
      'Hao, Y. et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nature Biotechnology 42, 293-304 (2024). https://doi.org/10.1038/s41587-023-01767-y',
      'Stuart, T. et al. Comprehensive Integration of Single-Cell Data. Cell 177, 1888-1902 (2019). https://doi.org/10.1016/j.cell.2019.05.031',
      'Korsunsky, I. et al. Fast, sensitive and accurate integration of single-cell data with Harmony. Nature Methods 16, 1289-1296 (2019). https://doi.org/10.1038/s41592-019-0619-0',
      'Chen, E.Y. et al. Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool. BMC Bioinformatics 14, 128 (2013). https://doi.org/10.1186/1471-2105-14-128',
      'Kuleshov, M.V. et al. Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. Nucleic Acids Research 44(W1), W90-W97 (2016). https://doi.org/10.1093/nar/gkw377',
      'Badia-i-Mompel, P. et al. decoupleR: ensemble of computational methods to infer biological activities from omics data. Bioinformatics Advances 2(1), vbac016 (2022). https://doi.org/10.1093/bioadv/vbac016',
      'Garcia-Alonso, L., Holland, C.H., Ibrahim, M.M., Turei, D. & Saez-Rodriguez, J. Benchmark and integration of resources for the estimation of human transcription factor activities. Genome Research 29, 1363-1375 (2019). https://doi.org/10.1101/gr.240663.118',
      'Schubert, M. et al. Perturbation-response genes reveal signaling footprints in cancer gene expression. Nature Communications 9, 20 (2018). https://doi.org/10.1038/s41467-017-02391-6'
    ]
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
    The report of merged single cell samples/plates<br/>
    EOS
    @required_columns = ['Name', 'Species', 'refBuild', 'refFeatureFile', 'Static Report']
    @required_params = []
    # optional params
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '60'
    @params['ram','description'] = 'use at least 20G per sample'
    @params['scratch'] = '100'
    @params['scratch', "context"] = "slurm"
    @params['node'] = ''
    @params['name'] = 'SCReportMultipleSamplesSeurat'
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "referfence genome assembly"
    @params['refFeatureFile'] = 'genes.gtf'
    @params['refFeatureFile', "context"] = "ScSeuratCombine"
    # --- Integration ---
    @params['integrationMethod', 'hr-header'] = "Integration"
    @params['integrationMethod'] = ['Harmony', 'CCA', 'RPCA', 'none']
    @params['integrationMethod', 'description'] = 'Harmony is the best general-purpose technique; use CCA for legacy reasons, use RPCA if the number of matching cells/cell types across your samples is small'
    @params['additionalFactors'] = ''
    @params['additionalFactors', 'description'] = "A comma-separated list of additional column names from the input dataset to use to label cells from a give sample. Useful for adding additional variables beyond 'Condition' and 'Batch' to the object. This information is also used by Harmony if Harmony is selected as the integration method. Use only the column name without '[Factor]'. Example: Patient,Tissue"
    @params['excludeGenePattern'] = ''
    @params['excludeGenePattern', 'description'] = "Regex of gene symbols to drop from the RNA assay of every sample before SCTransform, so they cannot be selected as variable features or drive the integration. Leave empty for normal runs. For a TCR/BCR-independent view of a lymphocyte dataset use: ^TR[ABGD][VDJC][0-9P]|^TR[ABGD]C[0-9]?$|^IG[HKL][VDJC]|^IGH[AEGMD][0-9]?$"
    # --- Normalization & Clustering ---
    @params['SCT.regress.CellCycle', 'hr-header'] = "Normalization & Clustering"
    @params['SCT.regress.CellCycle'] = false
    @params['SCT.regress.CellCycle', 'description'] = 'Choose CellCycle to be regressed out when using the SCTransform method if it is a bias.'
    @params['npcs'] = '30'
    @params['npcs', 'description'] = 'The maximal top dimensions (pcs) to use for reduction. Do not use more principal components than pcGenes (when used).'
    @params['pcGenes'] = ''
    @params['pcGenes', 'description'] = 'The genes used in supervised clustering'
    @params['resolution'] = '0.6'
    @params['resolution', 'description'] = 'Clustering resolution. A higher number will lead to more clusters.'
    # --- Cluster Markers ---
    @params['DE.method', 'hr-header'] = "Cluster Markers"
    @params['DE.method'] = ['wilcox', 'LR']
    @params['DE.method', 'description'] = "Method to be used when calculating gene cluster markers and differentially expressed genes between conditions. Use LR to take into account the Batch and/or CellCycle."
    @params['DE.regress'] = ['Batch', 'CellCycle']
    @params['DE.regress','multi_selection'] = true
    @params['DE.regress', 'description'] = "Variables to regress when calculating gene cluster markers and differentially expressed genes. Only used with the LR method."
    @params['min.pct'] = 0.1
    @params['min.pct', 'description'] = 'Used in calculating cluster markers: The minimum fraction of cells in either of the two tested populations.'
    @params['logfc.threshold'] = 0.25
    @params['logfc.threshold', 'description'] = 'Used in calculating cluster markers: Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells.'
    # --- Cell Type Annotation ---
    @params['tissue', 'hr-header'] = "Cell Type Annotation"
    @params['tissue'] = []
    @params['tissue','multi_selection'] = true
    @params['tissue','all_selected'] = true
    @params['tissue', 'multi_selection_size'] = 10
    tissue = {}
    CSV.foreach("/srv/GT/databases/scGeneSets/CellMarker_2.0-2023-09-27/Cell_marker_All_tissueList.txt", headers: true, col_sep: "\t") do |e|
      tissue[e["tissue_class"]] = true
    end
    @params['tissue'] = tissue.keys.sort
    @params['tissue', 'description'] = 'Select the tissues from the CellMarker2 database to identify celltypes using AUCell'
    @params['enrichrDatabase'] = ['Human_Gene_Atlas', 'Tabula_Sapiens', 'Azimuth_2023', 'PanglaoDB_Augmented_2021',
                                  'CellMarker_2024', 'HuBMAP_ASCTplusB_augmented_2022', 'Allen_Brain_Atlas_10x_scRNA_2021', 'Mouse_Gene_Atlas', 'Tabula_Muris', ]
    @params['enrichrDatabase','multi_selection'] = true
    @params['enrichrDatabase','all_selected'] = true
    @params['SingleR'] = ['none', 'BlueprintEncodeData (human)', 'DatabaseImmuneCellExpressionData (human)', 'HumanPrimaryCellAtlasData (human)',
                          'MonacoImmuneData (human)', 'NovershternHematopoieticData (human)', 'ImmGenData (mouse)', 'MouseRNAseqData (mouse)']
    @params['SingleR', 'description'] = "Use reference datasets from the celldex package to find marker-based celltype annotation with SingleR"
    # --- Additional Options ---
    @params['computePathwayTFActivity', 'hr-header'] = "Additional Options"
    @params['computePathwayTFActivity'] = false
    @params['computePathwayTFActivity', 'description'] = 'Whether to calculate the TF and pathway activities (Note: Only for human and mouse)'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Dev/R"]
    @inherit_columns = ["Order Id"]
  end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, '00index.html')
    seurat_file = File.join(report_file, 'scData.qs2')
    {'Name'=>@params['name'],
     'Species'=>(dataset = @dataset.first and dataset['Species']),
     'Static Report [Link]'=>report_link,
     'Report [File]'=>report_file,
     'SeuratObject [Link]'=>seurat_file
    }.merge(extract_columns(colnames: @inherit_columns))
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('refFeatureFile')
      @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    end
  end
  def commands
    #command = "module load #{@params["Rversion"]}\n"
    run_RApp("EzAppScSeuratCombine")
  end
end

if __FILE__ == $0

end
