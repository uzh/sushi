#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class ScSeuratApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'ScSeurat'
    # Seurat (v5, matches installed 5.5.1) unconditional. scDblFinder unconditional
    # (tryCatch-wrapped; only doublet removal gated on keepDoublets). emptyDrops
    # runs when raw matrix has extra barcodes (data-driven). cyclone/decontX/SoupX/
    # enrichR/AUCell/SingleR/celldex/decoupleR/dorothea/progeny/sc-type/
    # mLLMCelltype/CyteTypeR/AzimuthPanHuman all gated on their own params (species
    # Human/Mouse, estimateAmbient, enrichrDatabase, tissue, SingleR,
    # computePathwayTFActivity, sctype.enabled, mLLMCelltype, CyteTypeR,
    # AzimuthPanHuman respectively) -- listed regardless of gating.
    # NOTE: Azimuth::RunAzimuth, cellxgene_annotation, STACAS and schard exist in
    # the shared seuratUtils.R helpers but are unreachable here -- ScSeuratApp
    # never declares param$Azimuth/cellxgeneUrl/cellxgeneLabel/featSelectionMethod.
    @citation = [
      'Hao, Y. et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nature Biotechnology 42, 293-304 (2024). https://doi.org/10.1038/s41587-023-01767-y',
      'Germain, P.-L. et al. Doublet identification in single-cell sequencing data using scDblFinder. F1000Research 10, 979 (2022). https://doi.org/10.12688/f1000research.73600.2',
      'Lun, A.T.L. et al. EmptyDrops: distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data. Genome Biology 20, 63 (2019). https://doi.org/10.1186/s13059-019-1662-y',
      'Scialdone, A. et al. Computational assignment of cell-cycle stage from single-cell transcriptome data. Methods 85, 54-61 (2015). https://doi.org/10.1016/j.ymeth.2015.06.021',
      'Yang, S. et al. Decontamination of ambient RNA in single-cell RNA-seq with DecontX. Genome Biology 21, 57 (2020). https://doi.org/10.1186/s13059-020-1950-6',
      'Young, M.D. & Behjati, S. SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data. GigaScience 9(12), giaa151 (2020). https://doi.org/10.1093/gigascience/giaa151',
      'Chen, E.Y. et al. Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool. BMC Bioinformatics 14, 128 (2013). https://doi.org/10.1186/1471-2105-14-128',
      'Kuleshov, M.V. et al. Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. Nucleic Acids Research 44(W1), W90-W97 (2016). https://doi.org/10.1093/nar/gkw377',
      'Aibar, S. et al. SCENIC: single-cell regulatory network inference and clustering. Nature Methods 14, 1083-1086 (2017). https://doi.org/10.1038/nmeth.4463',
      'Aran, D. et al. Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. Nature Immunology 20(2), 163-172 (2019). https://doi.org/10.1038/s41590-018-0276-y',
      'Aran, D. et al. celldex: Reference Index for Cell Types. R package version 1.22.0. https://doi.org/10.18129/B9.bioc.celldex',
      'Badia-i-Mompel, P. et al. decoupleR: ensemble of computational methods to infer biological activities from omics data. Bioinformatics Advances 2(1), vbac016 (2022). https://doi.org/10.1093/bioadv/vbac016',
      'Garcia-Alonso, L., Holland, C.H., Ibrahim, M.M., Turei, D. & Saez-Rodriguez, J. Benchmark and integration of resources for the estimation of human transcription factor activities. Genome Research 29, 1363-1375 (2019). https://doi.org/10.1101/gr.240663.118',
      'Schubert, M. et al. Perturbation-response genes reveal signaling footprints in cancer gene expression. Nature Communications 9, 20 (2018). https://doi.org/10.1038/s41467-017-02391-6',
      'Ianevski, A., Giri, A.K. & Aittokallio, T. Fully-automated and ultra-fast cell-type identification using specific marker combinations from single-cell transcriptomic data. Nature Communications 13, 1246 (2022). https://doi.org/10.1038/s41467-022-28803-w',
      'Yang et al. Large language model consensus substantially improves the cell type annotation accuracy for scRNA-seq data. Communications Biology (2026). https://doi.org/10.1038/s42003-026-10420-8',
      'Ahuja, G. et al. Multi-agent AI enables evidence-based cell annotation in single-cell transcriptomics. bioRxiv (2025) [preprint, not peer-reviewed]. https://doi.org/10.1101/2025.11.06.686964',
      'Satija Lab. Pan-Human Azimuth. https://satijalab.org/pan_human_azimuth/ [preprint referenced on this page could not be independently verified]'
    ]
    @params['process_mode'] = 'SAMPLE'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
Single cell report<br/>
    EOS
    @required_columns = ['Name', 'Species', 'refBuild', 'CountMatrix', 'ResultDir', 'Condition']
    @required_params = ['name']
    # optional params
    @params['cores'] = '4'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '100'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '100'
    @params['scratch', "context"] = "slurm"
    @params['name'] = 'ScSeurat'
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "referfence genome assembly"
    @params['refFeatureFile'] = 'genes.gtf'
    @params['refFeatureFile', "context"] = "ScSeurat"
    @params['geneCountModel'] = ''
    @params['geneCountModel', 'description'] = '(STARsolo Input Only) The gene count model, i.e. Solo features, to use from the previous step'
    # --- Quality Control ---
    @params['nUMI', 'hr-header'] = "Quality Control"
    @params['nUMI'] = ''
    @params['nUMI', 'description'] = "Low quality cells have less than 'nUMI' UMIs. Only when applying fixed thresholds"
    @params['ngenes'] = ''
    @params['ngenes', 'description'] = "Low quality cells have less than 'ngenes' genes. Only when applying fixed thresholds"
    @params['perc_mito'] = ''
    @params['perc_mito', 'description'] = "Low quality cells have more than 'perc_mito' percent of mitochondrial genes. Only when applying fixed thresholds"
    @params['perc_riboprot'] = '70'
    @params['perc_riboprot', 'description'] = "Low quality cells have more than 'perc_riboprot' percent of ribosomal genes. Only when applying fixed thresholds"
    @params['cellsFraction'] = 0.0
    @params['cellsFraction', 'description'] = 'A gene will be kept if it is expressed in at least this fraction of cells'
    @params['geneMinUMI'] = 1
    @params['geneMinUMI', 'description'] = 'A gene will be kept if it has at least this many UMIs in the fraction of cells specified before'
    @params['filterByExpression'] = ''
    @params['filterByExpression', 'description'] = 'Keep cells according to specific gene expression. i.e. Set > 1 | Pkn3 > 1'
    @params['estimateAmbient'] = true
    @params['estimateAmbient', 'description'] = 'Run SoupX and DecontX to estimate ambient RNA levels'
    # --- Normalization & Clustering ---
    @params['SCT.regress.CellCycle', 'hr-header'] = "Normalization & Clustering"
    @params['SCT.regress.CellCycle'] = false
    @params['SCT.regress.CellCycle', 'description'] = 'Choose CellCycle to be regressed out when using the SCTransform method if it is a bias.'
    @params['npcs'] = 20
    @params['npcs', 'description'] = 'The maximal top dimensions (pcs) to use for reduction. Do not use more principal components than pcGenes (when used).'
    @params['pcGenes'] = ''
    @params['pcGenes', 'description'] = 'The genes used in supervised clustering'
    @params['resolution'] = [0.6, 0.2, 0.4, 0.6, 0.8, 1]
    @params['resolution', 'description'] = 'Clustering resolution. A higher number will lead to more clusters.'
    # --- Cluster Markers ---
    @params['DE.method', 'hr-header'] = "Cluster Markers"
    @params['DE.method'] = ['wilcox', 'LR']
    @params['DE.method', 'description'] ='Method to be used when calculating gene cluster markers. Use LR if you want to include cell cycle in the regression model.'
    @params['min.pct'] = 0.1
    @params['min.pct', 'description'] = 'Used in calculating cluster markers: The minimum fraction of cells in either of the two tested populations.'
    @params['min.diff.pct'] = 0.1
    @params['min.diff.pct', 'description'] = 'Used for filtering cluster markers: The minimum difference of cell fraction of the two tested populations.'
    @params['logfc.threshold'] = 0.25
    @params['logfc.threshold', 'description'] = 'Used in calculating cluster markers: Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells.'
    @params['pvalue_allMarkers'] = 0.01
    @params['pvalue_allMarkers', 'description'] = 'Used for filtering cluster markers: adjusted pValue threshold for marker detection.'
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
    @params['Azimuth'] = ["none", "adiposeref (human)", "bonemarrowref (human)", "fetusref (human)", "heartref (human)", "humancortexref (human)",
                          "kidneyref (human)", "lungref (human)", "pancreasref (human)", "pbmcref (human)", "tonsilref (human)", "/srv/GT/databases/Azimuth/humanLiver_Azimuth_v1.0 (human)", "/srv/GT/databases/Azimuth/murineLiver_v1 (mouse)",
                          "mousecortexref (mouse)"]
    @params['AzimuthPanHuman'] = false
    @params['AzimuthPanHuman', 'description'] = 'Enable Azimuth Pan-Human cell type annotation. If the upstream CellRanger run (>=10.1.0) already produced outs/cell_types/, those labels are reused and nothing leaves the site - verified bit-for-bit identical to the API result. Only when that file is absent does this fall back to the external CloudAzimuth API, which sends the expression matrix off-site. Human only - ezRun checks refBuild and skips non-human datasets automatically.'
    @params['SingleR'] = ['none', 'BlueprintEncodeData (human)', 'DatabaseImmuneCellExpressionData (human)', 'HumanPrimaryCellAtlasData (human)',
                          'MonacoImmuneData (human)', 'NovershternHematopoieticData (human)', 'ImmGenData (mouse)', 'MouseRNAseqData (mouse)']
    @params['SingleR', 'description'] = "Use reference datasets from the celldex package to find marker-based celltype annotation with SingleR"
    @params['sctype.enabled'] = true
    @params['sctype.enabled', 'description'] = 'Enable scType automatic cell type annotation (human and mouse supported)'
    @params['sctype.tissue'] = ["auto", "Immune system", "Liver", "Pancreas", "Kidney", "Eye", "Brain", "Lung", "Adrenal", "Heart", "Intestine", "Muscle", "Placenta", "Spleen", "Stomach", "Thymus"]
    @params['sctype.tissue', 'description'] = 'Tissue type for scType annotation. Select "auto" for automatic detection or specify the tissue type for more accurate results'
    @params['CyteTypeR'] = false
    @params['CyteTypeR', 'description'] = 'Enable CyteTypeR AI-powered cell type annotation (requires API key)'
    @params['CyteTypeR.apiKey'] = ''
    @params['CyteTypeR.apiKey', 'description'] = 'Nygen Analytics CyteType API key (or set CYTETYPE_API_KEY env var)'
    @params['CyteTypeR.studyContext'] = ''
    @params['CyteTypeR.studyContext', 'description'] = 'Optional biological context for better annotation (e.g., "Human PBMC from healthy donors"). Auto-generated if empty.'
    @params['mLLMCelltype'] = true
    @params['mLLMCelltype', 'description'] = 'Enable mLLMCelltype cluster annotation on the FGCZ-internal vLLM. Only the top marker gene NAMES per cluster are sent and the server is on-site, so no expression data leaves FGCZ and no API key is needed.'
    @params['mLLMCelltype.tissue'] = ''
    @params['mLLMCelltype.tissue', 'description'] = 'Tissue or sorted population, e.g. "PBMC", "sorted B cells", "lung", "CD8+ TILs". Strongly recommended: with no context the labels get much worse (on a sorted B-cell dataset, 4/14 clusters correct without it vs 13/14 with it).'
    @params['cellxgeneUrl'] = ''
    @params['cellxgeneUrl', 'description'] = 'Choose a download URL to a Seurat rds file of a dataset from here: https://cellxgene.cziscience.com/datasets'
    @params['cellxgeneLabel'] = ''
    @params['cellxgeneLabel', 'description'] = 'Metadata column for cell type labels. "cell_type" is a good starting choice'
    # --- Additional Options ---
    @params['computePathwayTFActivity', 'hr-header'] = "Additional Options"
    @params['computePathwayTFActivity'] = false
    @params['computePathwayTFActivity', 'description'] = 'Whether to calculate the TF and pathway activities (Note: Only for human and mouse)'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def preprocess
    @random_string = (1..12).map{[*('a'..'z')].sample}.join
  end
  def next_dataset
    report_file = File.join(@result_dir, "#{@dataset['Name']}_SCReport")
    report_link = File.join(report_file, '00index.html')
    # Forward CountMatrix + ResultDir so downstream apps (ScMultiOmics) can
    # locate sibling modalities (vdj_t/, atac_fragments.tsv.gz, Antibody
    # Capture in the H5) — these live next to the original CellRanger output,
    # not next to scData.qs2. Both columns are required inputs to ScSeurat
    # (see @required_columns line 18), so the values are guaranteed present.
    # Use [Link] (not [File]) for both: we are propagating REFERENCES to the
    # upstream CellRanger output, not producing new files. [File] would make
    # SUSHI's job_footer attempt `g-req copy <upstream_dir> <upstream_gstore>`,
    # which fails with "Destination path already exists".
    {'Name'=>@dataset['Name'],
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'refFeatureFile'=>@params['refFeatureFile'],
     'CountMatrix [Link]'=>@dataset['CountMatrix'],
     'ResultDir [Link]'=>@dataset['ResultDir'],
     'Static Report [Link]'=>report_link,
     'SC Cluster Report [File]'=>report_file,
     'SC Seurat [Link]'=>File.join(report_file, "scData.qs2"),
    }.merge(extract_columns(@inherit_tags))
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('refFeatureFile')
      @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    end
    if dataset_has_column?('soloFeatures')
      @params['geneCountModel'] = @dataset[0]['soloFeatures'].split(',')
    else
      @params.delete('geneCountModel')
    end
  end
  def commands
    #command = "module load #{@params["Rversion"]}\n"
    run_RApp("EzAppScSeurat")
  end
end

if __FILE__ == $0

end

