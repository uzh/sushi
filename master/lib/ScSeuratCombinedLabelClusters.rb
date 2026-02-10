#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class ScSeuratCombinedLabelClusters < SushiFabric::SushiApp
  def initialize
    super
    @name = 'LabelCombinedClusters'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
    The report of merged single cell samples/plates<br/>
    EOS
    @required_columns = ['Name', 'Species', 'Static Report', 'Report', 'SeuratObject']
    @required_params = ['ClusterAnnotationFile']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '60'
    @params['scratch'] = '100'
    @params['node'] = ''
    @params['name'] = 'SCReportMultipleSamplesSeurat'
    # --- Cluster Labels ---
    @params['ClusterAnnotationFile', 'hr-header'] = "Cluster Labels"
    @params['ClusterAnnotationFile'] = ''
    @params['ClusterAnnotationFile', 'file_upload'] = true
    @params['ClusterAnnotationFile', 'description'] = "A 3-column mapping of old cluster to new cluster labels in .xlsx format. Use the 'clusterInfos.xlsx' file as a template. The first column indicates the sample name. The second column are the old cluster labels. The third column are the new cluster labels. The first row should be a header indicating the column names. We recommend using the names provided in 'clusterInfos.xlsx', namely 'Sample', 'Cluster', and 'ClusterLabel'."
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
     'SeuratObject [Link]'=>seurat_file,
    }.merge(extract_columns(colnames: @inherit_columns))
  end
  def commands
    run_RApp("EzAppScSeuratCombinedLabelClusters")
  end
end

if __FILE__ == $0

end
