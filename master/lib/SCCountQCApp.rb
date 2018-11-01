#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class SCCountQCApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'SCCountQC'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
Quality control for singel cell alignment and counts<br/>
    EOS
    @required_columns = ['Name', 'Species', 'refBuild', 'CountMatrix', 'BAM']
    @required_params = ['name']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '100'
    @params['scratch'] = '200'
    @params['name'] = 'SCCount_QC'
    @params['refBuild'] = ref_selector
    @params['paired'] = false
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['refFeatureFile'] = 'genes.gtf'
    @params['featureLevel'] = ['gene', 'isoform']
    @params['transcriptTypes'] = ''
    @params['transcriptTypes', 'multi_selection'] = true
    @params['transcriptTypes', 'selected'] = 0
    @params['minReadsPerCell'] = 1500
    @params['minReadsPerCell', 'description'] = 'Filter cells with less reads counted on genes'
    @params['minReadsPerGene'] = 3
    @params['minReadsPerGene', 'description'] = 'Minimal number of reads per gene to be expressed'
    @params['minGenesPerCell'] = 500
    @params['minGenesPerCell', 'description'] = 'Filter cells with less genes expressed'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Dev/R", "Dev/jdk", "Tools/Picard", "Tools/samtools"]
  end
  def next_dataset
    report_file = File.join(@result_dir, "#{@dataset['Name']}_SCCountQC")
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@dataset['Name'],
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'refFeatureFile'=>@params['refFeatureFile'],
     'Static Report [Link]'=>report_link,
     'seuratOnline [Link]'=>"#{SHINY_SEURAT_ONLINE}?data=#{report_file}/sce.rds",
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
    run_RApp("EzAppSCCountQC")
  end
end

if __FILE__ == $0

end

