#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class VisiumHDSeuratApp < SushiFabric::SushiApp
  def initialize
    super
    @employee = true
    @name = 'VisiumHDSeuratApp'
    @params['process_mode'] = 'SAMPLE'
    @analysis_category = 'Spatial'
    @description =<<-EOS
Single cell report<br/>
    EOS
    @required_columns = ['Name', 'Species', 'refBuild', 'SpaceRangerDir', 'BinnedOutput']
    @required_params = ['name']
    # optional params
    @params['cores'] = '4'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '100'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '50'
    @params['scratch', "context"] = "slurm"
    @params['name'] = 'VisiumHDSeurat'
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "referfence genome assembly"
    @params['refFeatureFile'] = 'genes.gtf'
    @params['refFeatureFile', "context"] = "SeuratVisiumHD"
    @params['binSize'] = ['binned_outputs/square_008um', 'binned_outputs/square_016um', 'segmented_outputs']
    @params['binSize','description'] = 'Standard binSizes are 8 and 16. Other binSize are only available if the parameter --custom-bin-size was used in SpaceRanger'
    @params['DE.method'] = ['wilcox', 'LR']
    @params['DE.method', 'description'] ='Method to be used when calculating gene cluster markers. Use LR if you want to include cell cycle in the regression model.'
    @params['npcs'] = 50
    @params['npcs', 'description'] = 'The maximal dimensions to use for reduction. Do not use more pcs than pcGenes (when used).'
    @params['min.pct'] = 0.1
    @params['min.pct', 'description'] = 'Used in calculating cluster markers: The minimum fraction of cells in either of the two tested populations.'
    @params['logfc.threshold'] = 0.25
    @params['logfc.threshold', 'description'] = 'Used in calculating cluster markers: Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells.'
    @params['pcGenes'] = ''
    @params['pcGenes', 'description'] = 'The genes used in supvervised clustering'
    @params['clusterResolution'] = [2, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.5, 3]
    @params['clusterResolution', 'description'] = 'Clustering resolution. A higher number will lead to more clusters.'
    @params['lambda'] = '0.8'
    @params['lambda', 'description'] = 'BANKSY spatial weighting (0-1). Higher values give more weight to spatial neighbors.'
    @params['Niche_resolution'] = '0.5'
    @params['Niche_resolution', 'description'] = 'Resolution for BANKSY spatial niche clustering (higher = more niches)'
    @params['numis'] = '20'
    @params['numis', 'description'] = 'minimum number of UMIs required per spot, lower the threshold then using 08um bins or smaller'
    @params['ngenes'] = ''
    @params['ngenes', 'description'] = 'Low quality spots have less than "ngenes" genes. Only when applying fixed thresholds'
    @params['perc_mito'] = ''
    @params['perc_mito', 'description'] = 'Low quality spots have more than "perc_mito" percent of mitochondrial genes. Only when applying fixed thresholds'
    @params['perc_ribo'] = ''
    @params['perc_ribo', 'description'] = 'Low quality spots have more than "perc_ribo" percent of ribosomal genes. Only when applying fixed thresholds'
    @params['nmad'] = 3
    @params['nmad', 'description'] = 'Median absolute deviation (MAD) from the median value of each metric across all bins'
    @params['lambda'] = 0.8
    @params['lambda', 'description'] = 'BANKSY lambda: spatial weighting parameter (0-1). Larger values (0.8) find spatial domains; smaller values (0.2) perform cell typing.'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @params['Rversion'] = ["Dev/R/4.5.0"]
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
     'Seurat Visium HD [File]'=>report_file,
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
    command << run_RApp("EzAppVisiumHDSeurat")
  end
end

if __FILE__ == $0

end

