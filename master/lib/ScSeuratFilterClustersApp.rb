#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class ScSeuratFilterClustersApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'ScSeuratFilter'
    @params['process_mode'] = 'SAMPLE'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
Single cell report<br/>
    EOS
    @required_columns = ['Name', 'Species', 'refBuild', 'SC Seurat']
    @required_params = ['name']
    # optional params
    @params['cores'] = '4'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '60'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '100'
    @params['scratch', "context"] = "slurm"
    @params['name'] = 'ScSeurat'
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "referfence genome assembly"
    @params['refFeatureFile'] = 'genes.gtf'
    @params['refFeatureFile', "context"] = "ScSeuratFilterClusters"
    # --- Cluster Filtering ---
    @params['ClusterAnnotationFile', 'hr-header'] = "Cluster Filtering"
    @params['ClusterAnnotationFile'] = ''
    @params['ClusterAnnotationFile', 'file_upload'] = true
    @params['ClusterAnnotationFile', 'description'] = "A 3-column mapping of old cluster to new cluster labels in .xlsx format. Use the 'clusterInfos.xlsx' file as a template. The first column indicates the sample name. The second column are the old cluster labels. The third column indicates if the cluster should be kept or not. The first row should be a header indicating the column names. We recommend using the names provided in 'clusterInfos.xlsx', namely 'Sample', 'Cluster', and 'ClusterLabel'. To remove a cluster, label it as REMOVE."
    @params['identityColumn'] = 'seurat_clusters'
    @params['identityColumn', 'description'] = "The name of the identity column to use for calculating the markers and annotate the cells after the subsetting. By default, the column 'seurat_clusters', the variable of the new clusterings, is used. Change this option if you want to e.g. re-use previous labels"
    # --- Normalization & Clustering ---
    @params['SCT.regress.CellCycle', 'hr-header'] = "Normalization & Clustering"
    @params['SCT.regress.CellCycle'] = false
    @params['SCT.regress.CellCycle', 'description'] = 'Choose CellCycle to be regressed out when using the SCTransform method if it is a bias.'
    @params['npcs'] = 15
    @params['npcs', 'description'] = 'The maximal top dimensions (pcs) to use for reduction. Do not use more principal components than pcGenes (when used).'
    @params['resolution'] = [0.6, 0.2, 0.4, 0.6, 0.8, 1]
    @params['resolution', 'description'] = 'Clustering resolution. A higher number will lead to more clusters.'
    # --- Cluster Markers ---
    @params['DE.method', 'hr-header'] = "Cluster Markers"
    @params['DE.method'] = ['wilcox', 'LR']
    @params['DE.method', 'description'] ='Method to be used when calculating gene cluster markers. Use LR if you want to include cell cycle in the regression model.'
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
    @params['Azimuth'] = ["none", "adiposeref (human)", "bonemarrowref (human)", "fetusref (human)", "heartref (human)", "humancortexref (human)",
                          "kidneyref (human)", "lungref (human)", "pancreasref (human)", "pbmcref (human)", "tonsilref (human)", "/srv/GT/databases/Azimuth/humanLiver_Azimuth_v1.0 (human)",
                          "mousecortexref (mouse)"]
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
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def next_dataset
    report_file = File.join(@result_dir, "#{@dataset['Name']}_SCReport")
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@dataset['Name'],
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'refFeatureFile'=>@params['refFeatureFile'],
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
  end
  def commands
    run_RApp("EzAppScSeuratFilterClusters")
  end
end

if __FILE__ == $0

end

