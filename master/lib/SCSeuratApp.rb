#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class SCSeuratApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'SCSeurat'
    @params['process_mode'] = 'SAMPLE'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
Single cell Seurat clustering report<br/>
    EOS
    @required_columns = ['Name', 'Species', 'refBuild', 'CountMatrix', 'ResultDir']
    @required_params = ['name']
    # optional params
    @params['cores'] = '1'
    @params['ram'] = '30'
    @params['scratch'] = '50'
    @params['name'] = 'SCSeurat'
    @params['refBuild'] = ref_selector
    @params['paired'] = false
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['refFeatureFile'] = 'genes.gtf'
    @params['featureLevel'] = ['gene', 'isoform']
    @params['transcriptTypes'] = ''
    @params['transcriptTypes', 'multi_selection'] = true
    @params['transcriptTypes', 'selected'] = 0
    @params['scProtocol'] = ['10X', 'smart-Seq2']
    @params['minCellsPerGene'] = 3
    @params['minCellsPerGene', 'description'] = 'Minimum number of cells per gene for creating Seurat object'
    @params['minGenesPerCell'] = 1000
    @params['minGenesPerCell', 'description'] = 'Minimal number of genes per cell for Seurat filtering'
    @params['maxGenesPerCell'] = 7000
    @params['maxGenesPerCell', 'description'] = 'Maximal number of genes per cell for Seurat filtering'
    @params['maxMitoPercent'] = 25
    @params['maxMitoPercent', 'description'] = 'Maximal percentage of mitochondrial reads per cell for Seurat filtering'
    @params['pcs'] = 50
    @params['pcs', 'description'] = 'The maximal dimensions to use for reduction.'
    @params['vars.to.regress'] = ['nFeature_RNA', 'nCount_RNA', 'percent.mt']
    @params['vars.to.regress', 'multi_selection'] = true
    @params['vars.to.regress', 'selected'] = ['nFeature_RNA', 'percent.mt']
    @params['resolution'] = 0.5
    @params['resolution', 'description'] = 'Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.'
    @params['markersToShow'] = 10
    @params['markersToShow', 'description'] = 'The markers to show in the heatmap of cluster marker genes.'
    @params['knownMarkers'] = ''
    @params['knownMarkers', 'description'] = 'The markers to check in format of DC-like=Lgals3,Napsa;B cells=Cd79a,Ly6d; . Single and double quote are not allowed.'
    @params['runPseudoTime'] = false
    @params['runPseudoTime', 'description'] = 'Run PseudoTime for single cell data?'
    @params['all2allMarkers'] = false
    @params['all2allMarkers', 'description'] = 'Run all against all cluster comparisons?'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @params['Rversion'] = ["Dev/R/3.6.0", "Dev/R/3.5.1", "Dev/R/3.5.0", "Dev/R/3.4.2", "Dev/R/3.4.0", "Dev/R/3.3.0"]
    @modules = ["Dev/R"]
  end
  def preprocess
    @random_string = (1..12).map{[*('a'..'z')].sample}.join
  end
  def next_dataset
    report_file = File.join(@result_dir, "#{@dataset['Name']}_SCSeurat")
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@dataset['Name'],
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'refFeatureFile'=>@params['refFeatureFile'],
     'Static Report [Link]'=>report_link,
     'Live Report [Link]'=>"#{SHINY_EXPLORE_SCSEURAT}?data=#{report_file}/SCSeurat-#{@random_string}.rds",
     'Report [File]'=>report_file,
     'ResultDir [Link]'=>@dataset['ResultDir'],
    }
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('refFeatureFile')
      @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    end
    @params['transcriptTypes'] = @dataset[0]['transcriptTypes'].to_s.split(',')
    if dataset_has_column?('paired')
      @params['paired'] = @dataset[0]['paired']
    end
    if dataset_has_column?('strandMode')
      @params['strandMode'] = @dataset[0]['strandMode']
    end
  end
  def commands
    run_RApp("EzAppSCSeurat")
  end
end

if __FILE__ == $0

end

