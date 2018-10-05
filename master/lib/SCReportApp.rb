#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class SCReportApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'SCReport'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
Single cell report<br/>
    EOS
    @required_columns = ['Name', 'Species', 'refBuild', 'CountMatrix', 'BAM']
    @required_params = ['name']
    # optional params
    @params['cores'] = '4'
    @params['ram'] = '8'
    @params['scratch'] = '10'
    @params['name'] = 'SCReport'
    @params['refBuild'] = ref_selector
    @params['paired'] = false
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['refFeatureFile'] = 'genes.gtf'
    @params['featureLevel'] = ['gene', 'isoform']
    @params['transcriptTypes'] = ''
    @params['transcriptTypes', 'multi_selection'] = true
    @params['transcriptTypes', 'selected'] = 0
    @params['min_genes'] = 500
    @params['min_genes', 'description'] = 'Minimal number of genes for Seurat filtering'
    @params['max_genes'] = 3000
    @params['max_genes', 'description'] = 'Maximal number of genes for Seurat filtering'
    @params['min_counts'] = 50000
    @params['min_counts', 'description'] = 'Minimal counts for Seurat filtering'
    @params['pcs'] = 10
    @params['pcs', 'description'] = 'The maximal dimensions to use for reduction'
    @params['markersToCheck'] = ''
    @params['markersToCheck', 'description'] = 'The markers to check in format of DC-like=Lgals3,Napsa;B cells=Cd79a,Ly6d; '
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Dev/R"]
  end
  def next_dataset
    report_file = File.join(@result_dir, "#{@dataset['Name']}_SCReport")
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@dataset['Name'],
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'refFeatureFile'=>@params['refFeatureFile'],
     'Static Report [Link]'=>report_link,
     'Live Report [Link]'=>"#{SHINY_EXPLORE_SC}?data=#{report_file}/SCReport-#{@random_string}.RData",
     'Report [File]'=>report_file,
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
    run_RApp("EzAppSCReport")
  end
end

if __FILE__ == $0

end

