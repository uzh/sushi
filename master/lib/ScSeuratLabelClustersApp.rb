#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class ScSeuratLabelClustersApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'ScSeurat'
    @params['process_mode'] = 'SAMPLE'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
Single cell report<br/>
    EOS
    @required_columns = ['Name', 'Species', 'refBuild', 'SC Seurat']
    @required_params = ['name']
    # optional params
    @params['cores'] = '4'
    @params['ram'] = '60'
    @params['scratch'] = '100'
    @params['name'] = 'ScSeurat'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['ClusterAnnotationFile'] = ''
    @params['ClusterAnnotationFile', 'file_upload'] = true
    @params['ClusterAnnotationFile', 'description'] = "A 3-column mapping of old cluster to new cluster labels in .xlsx format. Use the 'clusterInfos.xlsx' file as a template. The first column indicates the sample name. The second column are the old cluster labels. The third column are the new cluster labels. The first row should be a header indicating the column names. We recommend using the names provided in 'clusterInfos.xlsx', namely 'Sample', 'Cluster', and 'ClusterLabel'."
    @params['DE.method'] = ['wilcox', 'LR']
    @params['DE.method', 'description'] ='Method to be used when calculating gene cluster markers. Use LR if you want to include cell cycle in the regression model.'
    @params['tissue'] = []
    @params['tissue','multi_selection'] = true
    @params['tissue','all_selected'] = true
    @params['tissue', 'multi_selection_size'] = 10
    tissue = {}
    CSV.foreach("/srv/GT/databases/scGeneSets/CellMarker_2.0-2023-09-27/Cell_marker_All_tissueList.txt", headers: true, col_sep: "\t") do |e|
      tissue[e["tissue_class"]] = true
    end
    @params['tissue'] = tissue.keys.sort
    @params['tissue', 'description'] = 'Select the tissues of origin in the CellMarker2 database.'
    @params['enrichrDatabase'] = ['Tabula_Muris', 'Tabula_Sapiens', 'Azimuth_Cell_Types_2021', 'PanglaoDB_Augmented_2021', 
                                  'CellMarker_Augmented_2021', 'Allen_Brain_Atlas_10x_scRNA_2021', 'Human_Gene_Atlas', 'Mouse_Gene_Atlas', ]
    @params['enrichrDatabase','multi_selection'] = true
    @params['enrichrDatabase','all_selected'] = true
    @params['Azimuth'] = ["none", "adiposeref", "bonemarrowref", "fetusref", "heartref", "humancortexref", 
                          "kidneyref", "lungref", "mousecortexref", "pancreasref", "pbmcref", "tonsilref"]
    @params['SingleR'] = ['none', 'BlueprintEncodeData', 'DatabaseImmuneCellExpressionData', 'HumanPrimaryCellAtlasData', 
                          'ImmGenData', 'MonacoImmuneData', 'MouseRNAseqData', 'NovershternHematopoieticData']
    @params['computePathwayTFActivity'] = true
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
     'SC Seurat'=>File.join(report_file, "scData.rds"),
    }.merge(extract_columns(@inherit_tags))
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('refFeatureFile')
      @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    end
  end
  def commands
    run_RApp("EzAppScSeuratLabelClusters")
  end
end

if __FILE__ == $0

end

