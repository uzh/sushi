#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class ScSeuratApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'ScSeurat'
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
    @params['SCT.regress.CellCycle'] = false
    @params['SCT.regress.CellCycle', 'description'] = 'Choose CellCycle to be regressed out when using the SCTransform method if it is a bias.'
    @params['DE.method'] = ['wilcox', 'LR']
    @params['DE.method', 'description'] ='Method to be used when calculating gene cluster markers. Use LR if you want to include cell cycle in the regression model.'
    @params['min.pct'] = 0.1
    @params['min.pct', 'description'] = 'Used in calculating cluster markers: The minimum fraction of cells in either of the two tested populations.'
    @params['min.diff.pct'] = 0.1
    @params['min.diff.pct', 'description'] = 'Used for filtering cluster markers: The minimum difference of cell fraction of the two tested populations.'
    @params['pvalue_allMarkers'] = 0.01
    @params['pvalue_allMarkers', 'description'] = 'Used for filtering cluster markers: adjusted pValue threshold for marker detection.'
    @params['logfc.threshold'] = 0.25
    @params['logfc.threshold', 'description'] = 'Used in calculating cluster markers: Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells.'
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
    @params['AzimuthPanHuman'] = false
    @params['AzimuthPanHuman', 'description'] = 'Enable Azimuth Pan-Human neural network-based cell type annotation (HUMAN DATASETS ONLY)'
    @params['AzimuthPanHuman.confidence.threshold'] = 0.5
    @params['AzimuthPanHuman.confidence.threshold', 'description'] = 'Confidence threshold for Azimuth Pan-Human annotation (0.0-1.0)'
    @params['SingleR'] = ['none', 'BlueprintEncodeData (human)', 'DatabaseImmuneCellExpressionData (human)', 'HumanPrimaryCellAtlasData (human)', 
                          'MonacoImmuneData (human)', 'NovershternHematopoieticData (human)', 'ImmGenData (mouse)', 'MouseRNAseqData (mouse)']
    @params['SingleR', 'description'] = "Use reference datasets from the celldex package to find marker-based celltype annotation with SingleR"
    @params['cellxgeneUrl'] = ''
    @params['cellxgeneUrl', 'description'] = 'Choose an download URL to an Seurat rds file of a dataset from here: https://cellxgene.cziscience.com/datasets'
    @params['cellxgeneLabel'] = ''
    @params['cellxgeneLabel', 'description'] = 'Specify the attribute of the dataset that should serve as cell type label'
    @params['sctype.enabled'] = true
    @params['sctype.enabled', 'description'] = 'Enable scType automatic cell type annotation (human and mouse supported)'
    @params['sctype.tissue'] = ["auto", "Immune system", "Liver", "Pancreas", "Kidney", "Eye", "Brain", "Lung", "Adrenal", "Heart", "Intestine", "Muscle", "Placenta", "Spleen", "Stomach", "Thymus"]
    @params['sctype.tissue', 'description'] = 'Tissue type for scType annotation. Select "auto" for automatic detection or specify the tissue type for more accurate results'
    @params['sctype.confidence.threshold'] = 0.25
    @params['sctype.confidence.threshold', 'description'] = 'Confidence threshold for scType annotation'
    @params['npcs'] = 20
    @params['npcs', 'description'] = 'The maximal top dimensions (pcs) to use for reduction. Do not use more principal components than pcGenes (when used).'
    @params['pcGenes'] = ''
    @params['pcGenes', 'description'] = 'The genes used in supervised clustering'
    @params['resolution'] = [0.6, 0.2, 0.4, 0.6, 0.8, 1]
    @params['resolution', 'description'] = 'Clustering resolution. A higher number will lead to more clusters.'
    @params['nUMI'] = ''
    @params['nUMI', 'description'] = "Low quality cells have less than 'nUMI' UMIs. Only when applying fixed thresholds"
    @params['ngenes'] = ''
    @params['ngenes', 'description'] = "Low quality cells have less than 'ngenes' genes. Only when applying fixed thresholds"
    @params['perc_mito'] = ''
    @params['perc_mito', 'description'] = "Low quality cells have more than 'perc_mito' percent of mitochondrial genes. Only when applying fixed thresholds"
    @params['perc_riboprot'] = '70'
    @params['perc_riboprot', 'description'] = "Low quality cells have more than 'perc_ribo' percent of ribosomal genes. Only when applying fixed thresholds"
    @params['cellsFraction'] = 0.0
    @params['cellsFraction', 'description'] = 'A gene will be kept if it is expressed in at least this fraction of cells'
    @params['geneMinUMI'] = 1
    @params['geneMinUMI', 'description'] = 'A gene will be kept if it has at least this many UMIs in the fraction of cells specified before'
    @params['filterByExpression'] = ''
    @params['filterByExpression', 'description'] = 'Keep cells according to specific gene expression. i.e. Set > 1 | Pkn3 > 1'
    @params['estimateAmbient'] = true
    @params['estimateAmbient', 'description'] = 'run SoupX and DecontX to estimate ambientRNA levels'
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

