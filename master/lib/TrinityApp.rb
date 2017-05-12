#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class TrinityApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Trinity'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Assemble'
    @description =<<-EOS
    RNA-seq de novo Assembly<br/>
<a href='https://github.com/trinityrnaseq/trinityrnaseq/wiki'>Trinity</a><br/>
    EOS
    @required_columns = ['Name', 'Read1', 'Species']
    @required_params = ['name']
    # optional params
    @params['cores'] = '12'
    @params['ram'] = '220'
    @params['scratch'] = '1000'
    @params['name'] = "Trinity_Assembly"
    @params['paired'] = false
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['trimAdapter'] = true
    @params['trimLeft'] = 5 
    @params['trimRight'] = 5
    @params['minTailQuality'] = 20
    @params['minAvgQuality'] = 20
    @params['minReadLength'] = 36
    @params['trinityOpt'] = '--min_kmer_cov 2'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Assembly/Trinity", "QC/Trimmomatic", "QC/Flexbar", "Aligner/Salmon"]
  end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  end
  def set_default_parameters
    @params['paired'] = dataset_has_column?('Read2')
    if dataset_has_column?('strandMode')
      @params['strandMode'] = @dataset[0]['strandMode']
    end
  end
  def next_dataset
    {'Name'=>@params['name'],
     'Fasta [File]'=>File.join(@result_dir, "#{@params['name']}.fasta"),
	  'Stats [File]'=>File.join(@result_dir, "assembly_stats.txt"),
	  'Abundance Counts [File]'=>File.join(@result_dir, "abundance_counts.txt"),
	  'Abundance TPM [File]'=>File.join(@result_dir, "abundance_TPM.txt"),
	  'Abundance TMM [File]'=>File.join(@result_dir, "abundance_TMM.txt"),
	  'ExN50 [File]'=>File.join(@result_dir, "ExN50_stats.txt")
    }
  end
  def commands
    run_RApp('EzAppTrinity')
  end
end

if __FILE__ == $0
  #usecase = EdgeRApp.new

  usecase.project = "p1001"
  usecase.user = 'masamasa'

  # set user parameter
  # for GUI sushi
  #usecase.params['process_mode'].value = 'SAMPLE'
  usecase.params['refBuild'] = 'mm10'
  usecase.params['paired'] = true
  usecase.params['strandMode'] = 'both'
  usecase.params['cores'] = 8
  usecase.params['node'] = 'fgcz-c-048'

  # also possible to load a parameterset csv file
  # mainly for CUI sushi
  #usecase.parameterset_tsv_file = 'tophat_parameterset.tsv'
  #usecase.parameterset_tsv_file = 'test.tsv'

  # set input dataset
  # mainly for CUI sushi
  #usecase.dataset_tsv_file = 'tophat_dataset.tsv'

  # also possible to load a input dataset from Sushi DB
  usecase.dataset_sushi_id = 3

  # run (submit to workflow_manager)
  usecase.run
  #usecase.test_run

end

