#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class SeuratXeniumApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'SeuratXenium'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Spatial'
    @description =<<-EOS
Seurat Analysis for Xenium Spatial Transcriptomics.<br/>
Includes QC, Normalization, Clustering, and RCTD Annotation.
    EOS
    @required_columns = ['Name','XeniumPath']
    @required_params = ['name']
    @params['cores'] = '8'
    @params['ram'] = '200'
    @params['scratch'] = '200'
    @params['name'] = 'SeuratXenium'
    @params['minCounts'] = '10'
    @params['minCounts', 'description'] = 'Minimum counts per cell for QC filtering'
    @params['minFeatures'] = '5'
    @params['minFeatures', 'description'] = 'Minimum features per cell for QC filtering'
    @params['Cluster_resolution'] = '0.5'
    @params['Cluster_resolution', 'description'] = 'Resolution for Seurat clustering (higher = more clusters)'
    @params['lambda'] = '0.8'
    @params['lambda', 'description'] = 'BANKSY spatial weighting (0-1). Higher values give more weight to spatial neighbors.'
    @params['Niche_resolution'] = '0.5'
    @params['Niche_resolution', 'description'] = 'Resolution for BANKSY spatial niche clustering (higher = more niches)'
    @params['rctdReference'] = ['None',
      'allen (mouse)',
      'azimuth (mouse)',
      'tabula_muris_senis (mouse)',
      'archmap (human)',
      'celltypist (human)',
      'disco (human)'
    ]
    @params['rctdReference', 'description'] = 'RCTD Reference folder. WARNING: RCTD requires 200GB-1TB RAM depending on dataset size. Available .rds files will be listed in the report.'
    @params['rctdFile'] = ''
    @params['rctdFile', 'description'] = 'Specific .rds file within reference folder (leave empty to see available files)'
    @params['rctdUMImin'] = '100'
    @params['rctdUMImin', 'description'] = 'Minimum UMI count for RCTD annotation. Cells below this threshold will not be classified.'
    @modules = ["Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def next_dataset
    report_dir = File.join(@result_dir, @params['name'])
    {'Name'=>@params['name'],
     'ReportData [File]'=>report_dir,
     'Report [Link]'=>File.join(report_dir, '00index.html'),
     'Species'=>(dataset = @dataset.first and dataset['Species'])
    }
  end
  def commands
    run_RApp("EzAppSeuratXenium")
  end
end
