#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class SpatialSeuratApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'SpatialSeurat'
    @params['process_mode'] = 'SAMPLE'
    @analysis_category = 'Spatial'
    @description =<<-EOS
Single cell report<br/>
    EOS
    @required_columns = ['Name', 'Species', 'refBuild', 'CountMatrix', 'ResultDir']
    @required_params = ['name']
    # optional params
    @params['cores'] = '4'
    @params['ram'] = '20'
    @params['scratch'] = '50'
    @params['name'] = 'SpatialSeurat'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['spotClean'] = false
    @params['spotClean', 'description'] = "run R Package spotClean on the raw count output of SpaceRanger instead of using the filtered output of SpaceRanger"
    @params['SCT.regress.CellCycle'] = false
    @params['SCT.regress.CellCycle', 'description'] = 'Choose CellCycle to be regressed out when using the SCTransform method if it is a bias.'
    @params['DE.method'] = ['wilcox', 'LR']
    @params['DE.method', 'description'] ='Method to be used when calculating gene cluster markers. Use LR if you want to include cell cycle in the regression model.'
    @params['npcs'] = 20
    @params['npcs', 'description'] = 'The maximal dimensions to use for reduction. Do not use more pcs than pcGenes (when used).'
    @params['min.pct'] = 0.1
    @params['min.pct', 'description'] = 'Used in calculating cluster markers: The minimum fraction of cells in either of the two tested populations.'
    @params['logfc.threshold'] = 0.25
    @params['logfc.threshold', 'description'] = 'Used in calculating cluster markers: Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells.'
    @params['pcGenes'] = ''
    @params['pcGenes', 'description'] = 'The genes used in supvervised clustering'
    @params['resolution'] = [0.6, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 2.5, 3]
    @params['resolution', 'description'] = 'Clustering resolution. A higher number will lead to more clusters.'
    @params['nreads'] = ''
    @params['nreads', 'description'] = 'Low quality spots have less than "nreads" reads. Only when applying fixed thresholds'
    @params['ngenes'] = ''
    @params['ngenes', 'description'] = 'Low quality spots have less than "ngenes" genes. Only when applying fixed thresholds'
    @params['perc_mito'] = ''
    @params['perc_mito', 'description'] = 'Low quality spots have more than "perc_mito" percent of mitochondrial genes. Only when applying fixed thresholds'
    @params['perc_ribo'] = ''
    @params['perc_ribo', 'description'] = 'Low quality spots have more than "perc_ribo" percent of ribosomal genes. Only when applying fixed thresholds'
    @params['nmad'] = 3
    @params['nmad', 'description'] = 'Median absolute deviation (MAD) from the median value of each metric across all spots'
    @params['cellsFraction'] = 0.0
    @params['cellsFraction', 'description'] = 'A gene will be kept if it is expressed in at least this fraction of spots'
    @params['nUMIs'] = 1
    @params['nUMIs', 'description'] = 'A gene will be kept if it has at least nUMIs in the fraction of spots specified before'
    @params['enrichrDatabase'] = ['Tabula_Muris', 'Tabula_Sapiens', 'Azimuth_Cell_Types_2021', 'PanglaoDB_Augmented_2021',
                                  'CellMarker_Augmented_2021', 'Allen_Brain_Atlas_10x_scRNA_2021', 'Human_Gene_Atlas', 'Mouse_Gene_Atlas']
    @params['enrichrDatabase','multi_selection'] = true
    @params['enrichrDatabase','all_selected'] = true
    @params['Azimuth'] = ["none", "adiposeref", "bonemarrowref", "fetusref", "heartref", "humancortexref",
                          "kidneyref", "lungref", "mousecortexref", "pancreasref", "pbmcref", "tonsilref","/srv/GT/databases/Azimuth/humanLiver_Azimuth_v1.0"]
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @params['Rversion'] = ["Dev/R/4.5.0","Dev/R/4.4.2"]
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
     'Report [File]'=>report_file,
    }.merge(extract_columns(@inherit_tags))
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('refFeatureFile')
      @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    end
  end
  def commands
    command = "module load #{@params["Rversion"]}\n"
    command << run_RApp("EzAppSpatialSeurat")
  end
end

if __FILE__ == $0

end

