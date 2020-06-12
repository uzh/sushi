#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class SCScranApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'SCScran'
    @params['process_mode'] = 'SAMPLE'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
Single cell report based on Scran<br/>
    EOS
    @required_columns = ['Name', 'Species', 'refBuild', 'CountMatrix', 'ResultDir']
    @required_params = ['name']
    # optional params
    @params['cores'] = '1'
    @params['ram'] = '30'
    @params['scratch'] = '50'
    @params['name'] = 'SCScran'
    @params['refBuild'] = ref_selector
    @params['paired'] = false
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['refFeatureFile'] = 'genes.gtf'
    @params['featureLevel'] = ['gene', 'isoform']
    @params['transcriptTypes'] = ''
    @params['transcriptTypes', 'multi_selection'] = true
    @params['transcriptTypes', 'selected'] = 0
    @params['scProtocol'] = ['10X', 'Smart-seq2']
    @params['snnK'] = '15'
    @params['snnK', 'description'] = 'An integer scalar specifying the number of nearest neighbors to consider during graph construction. Larger value, fewer clusters. Used in 10X.'
    @params['visMethod'] = ['TSNE', 'UMAP', 'DiffusionMap']
    @params['visMethod', 'description'] = 'Select the low-dimensional visusalization method.'
    @params['knownMarkers'] = ''
    @params['knownMarkers', 'description'] = 'The markers to check in format of DC-like=Lgals3,Napsa;B cells=Cd79a,Ly6d; . Single and double quote are not allowed.'
    @params['runPseudoTime'] = false
    @params['runPseudoTime', 'description'] = 'Run PseudoTime for single cell data?'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Dev/R"]
  end
  def preprocess
    @random_string = (1..12).map{[*('a'..'z')].sample}.join
  end
  def next_dataset
    report_file = File.join(@result_dir, "#{@dataset['Name']}_SCScran")
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@dataset['Name'],
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'refFeatureFile'=>@params['refFeatureFile'],
     'Static Report [Link]'=>report_link,
     'Live Report [Link]'=>"#{SHINY_EXPLORE_Scran}?data=#{report_file}/sce-#{@random_string}.rds",
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
    run_RApp("EzAppSCScran")
  end
end

if __FILE__ == $0

end

