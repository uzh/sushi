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
    @analysis_category = 'SpatialTrx'
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
    @params['SCT.regress'] = ['none', 'CellCycle']
    @params['SCT.regress', 'description'] = 'Choose CellCycle to be regressed out when using the SCTransform method if it is a bias.'
    @params['DE.method'] = ['wilcox', 'LR']
    @params['DE.method', 'description'] ='Method to be used when calculating gene cluster markers. Use LR if you want to include cell cycle in the regression model.'
    @params['npcs'] = 20
    @params['npcs', 'description'] = 'The maximal dimensions to use for reduction. Do not use more pcs than pcGenes (when used).'
    @params['pcGenes'] = ''
    @params['pcGenes', 'description'] = 'The genes used in supvervised clustering'
    @params['resolution'] = [0.6, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
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
    @params['cellsFraction'] = 0.05
    @params['cellsFraction', 'description'] = 'A gene will be kept if it is expressed in at least this fraction of spots'
    @params['nUMIs'] = 1
    @params['nUMIs', 'description'] = 'A gene will be kept if it has at least nUMIs in the fraction of spots specified before'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @params['Rversion'] = ["Dev/R/4.2.2", "Dev/R/4.2.0", "Dev/R/4.1.2", "Dev/R/4.1.0", "Dev/R/4.0.4", "Dev/R/4.0.3"]
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
    }
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

