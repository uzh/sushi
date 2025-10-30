#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class RnaBamStatsApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'RNA BamStats'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'QC'
    @description =<<-EOS
Quality control after the alignment of RNAseq reads<br/>
    EOS
    @required_columns = ['Name','BAM','BAI', 'refBuild', 'Species']
    @required_params = ['name', 'paired']
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '100'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '100'
    @params['scratch', "context"] = "slurm"
    @params['paired'] = false
    @params['paired', "context"] = "RnaBamStats"
    @params['name'] = 'RNA_BAM_Statistics'
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "referfence genome assembly"
    @params['refFeatureFile'] = 'genes.gtf'
    @params['refFeatureFile', "context"] = "RnaBamStats"
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['strandMode', "context"] = "RnaBamStats"
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Tools/samtools", "Dev/Python", "Tools/Picard", "Tools/BamUtil", "Dev/R"]
    @inherit_columns = ["Order Id"]
  end
  def next_dataset
    report_dir = File.join(@result_dir, @params['name'])
    {'Name'=>@params['name'],
     'Report [File]'=>report_dir,
     'Html [Link]'=>File.join(report_dir, '00index.html'),
     'Species'=>(dataset = @dataset.first and dataset['Species']),
     'refBuild'=>@params['refBuild'],
     'refFeatureFile'=>@params['refFeatureFile']
    }.merge(extract_columns(colnames: @inherit_columns))
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('refFeatureFile')
      @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    end
    if dataset_has_column?('paired')
      @params['paired'] = @dataset[0]['paired']
    end
    if dataset_has_column?('strandMode')
      @params['strandMode'] = @dataset[0]['strandMode']
    end
  end

  def commands
    run_RApp("EzAppRnaBamStats")
  end
end

if __FILE__ == $0
  usecase = BamStatsApp.new

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
