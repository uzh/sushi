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
    @params['ram'] = '60'
    @params['scratch'] = '100'
    @params['name'] = 'ScSeurat'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['geneCountModel'] = ''
    @params['geneCountModel', 'description'] = '(STARsolo Input Only) The gene count model, i.e. Solo features, to use from the previous step'
    @params['SCT.regress.CellCycle'] = false
    @params['SCT.regress', 'description'] = 'Choose CellCycle to be regressed out when using the SCTransform method if it is a bias.'
    @params['DE.method'] = ['wilcox', 'LR']
    @params['DE.method', 'description'] ='Method to be used when calculating gene cluster markers. Use LR if you want to include cell cycle in the regression model.'
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
    @params['npcs'] = 20
    @params['npcs', 'description'] = 'The maximal top dimensions (pcs) to use for reduction. Do not use more principal components than pcGenes (when used).'
    @params['pcGenes'] = ''
    @params['pcGenes', 'description'] = 'The genes used in supervised clustering'
    @params['resolution'] = [0.6, 0.2, 0.4, 0.6, 0.8, 1]
    @params['resolution', 'description'] = 'Clustering resolution. A higher number will lead to more clusters.'
    @params['nreads'] = ''
    @params['nreads', 'description'] = 'Low quality cells have less than "nreads" reads. Only when applying fixed thresholds'
    @params['ngenes'] = ''
    @params['ngenes', 'description'] = 'Low quality cells have less than "ngenes" genes. Only when applying fixed thresholds'
    @params['perc_mito'] = ''
    @params['perc_mito', 'description'] = 'Low quality cells have more than "perc_mito" percent of mitochondrial genes. Only when applying fixed thresholds'
    @params['perc_riboprot'] = '70'
    @params['perc_riboprot', 'description'] = 'Low quality cells have more than "perc_ribo" percent of ribosomal genes. Only when applying fixed thresholds'
    @params['cellsFraction'] = 0.0001
    @params['cellsFraction', 'description'] = 'A gene will be kept if it is expressed in at least this fraction of cells'
    @params['nUMIs'] = 1
    @params['nUMIs', 'description'] = 'A gene will be kept if it has at least nUMIs in the fraction of cells specified before'
    @params['filterByExpression'] = ''
    @params['filterByExpression', 'description'] = 'Keep cells according to specific gene expression. i.e. Set > 1 | Pkn3 > 1'
    @params['estimateAmbient'] = true
    @params['estimateAmbient', 'description'] = 'run SoupX and DecontX to estimate ambientRNA levels'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    #@params['Rversion'] = ["Dev/R/4.1.2", "Dev/R/4.1.0", "Dev/R/4.0.4", "Dev/R/4.0.3"]
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
     'SC Seurat'=>File.join(report_file, "scData.rds"),
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

