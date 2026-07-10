#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class DnaBamStatsApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'DNA BamStats'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'QC'
    @description = <<-EOS
Quality control after the alignment of DNA reads with one dataset-level report across all samples.<br/>
    EOS
    @required_columns = ['Name', 'BAM', 'BAI', 'refBuild', 'Species', 'Read Count']
    @required_params = ['name', 'paired']
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '50'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '100'
    @params['scratch', "context"] = "slurm"
    @params['paired'] = false
    @params['paired', "context"] = "DnaBamStats"
    @params['name'] = 'DNA_BAM_Statistics'
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "referfence genome assembly"
    @params['refFeatureFile'] = 'genes.gtf'
    @params['refFeatureFile', "context"] = "DnaBamStats"
    @params['pixelDist'] = '2500'
    @params['pixelDist', "context"] = "DnaBamStats"
    @params['runQualimap'] = true
    @params['runQualimap', "context"] = "DnaBamStats"
    @params['runPicard'] = true
    @params['runPicard', "context"] = "DnaBamStats"
    @params['specialOptions'] = ''
    @params['mail'] = ""
    # No explicit jdk here: QC/Qualimap loads Dev/jdk/8 and Tools/Picard loads
    # Dev/jdk/21, so Picard (listed after Qualimap) leaves jdk21 as the ambient
    # java. Picard needs jdk21; Qualimap tolerates it because the R app sets
    # JAVA_OPTS to skip its Java-8-only -XX:MaxPermSize flag. Keep QC/Qualimap
    # before Tools/Picard so the ambient java stays jdk21.
    @modules = ["Tools/samtools", "QC/Qualimap", "Tools/Picard", "Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end

  def next_dataset
    report_dir = File.join(@result_dir, @params['name'])
    {
      'Name' => @params['name'],
      'Report [File]' => report_dir,
      'Html [Link]' => File.join(report_dir, '00index.html'),
      'Species' => (dataset = @dataset.first and dataset['Species']),
      'refBuild' => @params['refBuild'],
      'refFeatureFile' => @params['refFeatureFile']
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
    run_RApp("EzAppDnaBamStats")
  end
end

if __FILE__ == $0
  usecase = DnaBamStatsApp.new

  usecase.project = "p1001"
  usecase.user = "masa"

  # set user parameter
  # for GUI sushi
  # usecase.params['process_mode'].value = 'DATASET'
  # usecase.params['refBuild'] = 'TAIR10'
  # usecase.params['paired'] = true
  # usecase.params['cores'] = 2
  # usecase.params['node'] = 'fgcz-c-048'

  # also possible to load a parameterset csv file
  # mainly for CUI sushi
  usecase.parameterset_tsv_file = 'tophat_parameterset.tsv'

  # set input dataset
  # mainly for CUI sushi
  usecase.dataset_tsv_file = 'tophat_dataset.tsv'

  # also possible to load a input dataset from Sushi DB
  # usecase.dataset_sushi_id = 1

  # run (submit to workflow_manager)
  usecase.run
  # usecase.test_run
end
