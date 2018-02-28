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
    @required_columns = ['Name', 'Species', 'featureLevel', 'refFeatureFile', 'CellDataset', 'CountMatrix', 'BAM', 'STARLog']
    @required_params = ['refBuild']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '100'
    @params['name'] = 'SCCount_QC'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['featureLevel'] = ['gene', 'isoform']
    @params['transcriptTypes'] = ''
    @params['transcriptTypes', 'multi_selection'] = true
    @params['transcriptTypes', 'selected'] = 0
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
     'Static Report [Link]'=>report_link,
     'Report [File]'=>report_file,
    }
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('refFeatureFile')
      @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    end
    @params['transcriptTypes'] = @dataset[0]['transcriptTypes'].to_s.split(',')
  end
  def commands
    run_RApp("EzAppSCCountQC")
  end
end

if __FILE__ == $0

end

