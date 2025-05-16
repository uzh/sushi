#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class ScSeuratCombineApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'ScSeuratCombine'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
    The report of merged single cell samples/plates<br/>
    EOS
    @required_columns = ['Name', 'Species', 'refBuild', 'refFeatureFile', 'Static Report']
    @required_params = []
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '60'
    @params['ram','description'] = 'use at least 20G per sample'
    @params['scratch'] = '100'
    @params['node'] = ''
    @params['name'] = 'SCReportMultipleSamplesSeurat'
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
    @params['enrichrDatabase'] = ['Tabula_Muris', 'Tabula_Sapiens', 'Azimuth_Cell_Types_2021', 'PanglaoDB_Augmented_2021', 'CellMarker_Augmented_2021', 'Allen_Brain_Atlas_10x_scRNA_2021', 'Human_Gene_Atlas', 'Mouse_Gene_Atlas', ]
    @params['enrichrDatabase','multi_selection'] = true
    @params['enrichrDatabase','all_selected'] = true
    @params['SingleR'] = ['none', 'BlueprintEncodeData', 'DatabaseImmuneCellExpressionData', 'HumanPrimaryCellAtlasData', 
                          'ImmGenData', 'MonacoImmuneData', 'MouseRNAseqData', 'NovershternHematopoieticData']
    @params['pcGenes'] = ''
    @params['pcGenes', 'description'] = 'The genes used in supvervised clustering'
    @params['npcs'] = '30'
    @params['npcs', 'description'] = 'Number of principal components to use for dimensionality reduction. Do not use more pcs than pcGenes (when used).'
    @params['resolution'] = '0.6'
    @params['resolution', 'description'] = 'Value between 0 and 1. A higher value will lead to larger communities.'
    @params['SCT.regress.CellCycle'] = false
    @params['SCT.regress.CellCycle', 'description'] = "Variable to regress when processing the counts with the SCTransform method."
    @params['DE.method'] = ['wilcox', 'LR']
    @params['DE.method', 'description'] = "Method to be used when calculating gene cluster markers and differentially expressed genes between conditions."
    @params['DE.regress'] = ['Batch', 'CellCycle']
    @params['DE.regress','multi_selection'] = true
    @params['DE.regress', 'description'] = "Variables to regress when calculating gene cluster markers and differentially expressed genes. Only used with the LR method."
    @params['min.pct'] = 0.1
    @params['min.pct', 'description'] = 'Used in calculating cluster markers: The minimum fraction of cells in either of the two tested populations.'
    @params['logfc.threshold'] = 0.25
    @params['logfc.threshold', 'description'] = 'Used in calculating cluster markers: Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells.'
    @params['integrationMethod'] = ['CCA', 'RPCA', 'Harmony', 'STACAS', 'none']
    @params['integrationMethod', 'description'] = 'use RPCA if the number of matching cells/cell types across your samples is small, STACAS if you have prior information on the labels'
    @params['additionalFactors'] = ''
    @params['additionalFactors', 'description'] = "A comma-separated list of additional column names from the input dataset to use to label cells from a give sample. Useful for adding additional variables beyond 'Condition' and 'Batch' to the object. This information is also used by Harmony if Harmony is selected as the integration method. Use only the column name without '[Factor]'. Example: Patient,Tissue"
    @params['STACASAnnotationFile'] = ''
    @params['STACASAnnotationFile', 'file_upload'] = true
    @params['STACASAnnotationFile', 'description'] = "(STACAS integration only) A 3-column mapping of old cluster to new cluster labels in .xlsx format. Use the 'clusterInfos.xlsx' file as a template. The first column indicates the sample name. The second column are the old cluster labels. The third column are the new cluster labels. The first row should be a header indicating the column names. We recommend using the names provided in 'clusterInfos.xlsx', namely 'Sample', 'Cluster', and 'ClusterLabel'. If you do not know the identity of a cluster with certainty, mark it 'unknown', rather than giving it your best guess. STACAS is designed to handle partial annotations."
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
