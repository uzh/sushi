#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class DnaQCApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'DNA QC'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'QC'
    @description =<<-EOS
Quality control after the alignment of RNAseq reads<br/>
    EOS
    @required_columns = ['Name','BAM','BAI', 'refBuild', 'Species']
    @required_params = ['name', 'paired']
    @params['cores'] = '8'
    @params['ram'] = '50'
    @params['scratch'] = '100'
    @params['paired'] = false
    @params['name'] = 'DNA_QC_Statistics'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Tools/samtools", "QC/Qualimap", "Tools/Picard", "Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def next_dataset
    report_dir = File.join(@result_dir, @params['name'])
    {'Name'=>@params['name'],
     'Report [File]'=>report_dir,
     'Html [Link]'=>File.join(report_dir, '00index.html'),
     'Species'=>(dataset = @dataset.first and dataset['Species']),
     'refBuild'=>@params['refBuild'],
     'refFeatureFile'=>@params['refFeatureFile']
    }.merge(extract_columns(@inherit_tags))
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('refFeatureFile')
      @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    end
    if dataset_has_column?('paired')
      @params['paired'] = @dataset[0]['paired']
    end
  end

  def commands
    run_RApp("EzAppDnaQC")
  end
end

if __FILE__ == $0
  usecase = DnaQCApp.new

  usecase.project = "p1001"
  usecase.user = "masa"

  # set user parameter
  # for GUI sushi
  #usecase.params['process_mode'].value = 'SAMPLE'
  #usecase.params['refBuild'] = 'TAIR10'
  #usecase.params['paired'] = true
  #usecase.params['cores'] = 2
  #usecase.params['node'] = 'fgcz-c-048'

  # also possible to load a parameterset csv file
  # mainly for CUI sushi
  usecase.parameterset_tsv_file = 'tophat_parameterset.tsv'
  #usecase.params['name'] = 'name'

  # set input dataset
  # mainly for CUI sushi
  usecase.dataset_tsv_file = 'tophat_dataset.tsv'

  # also possible to load a input dataset from Sushi DB
  #usecase.dataset_sushi_id = 1

  # run (submit to workflow_manager)
  usecase.run
  #usecase.test_run

end
