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
    @required_columns = ['Name','Read1', 'Species', ]
    @required_params = ['name']
    # optional params
    @params['cores'] = '12'
    @params['ram'] = '220'
    @params['scratch'] = '1000'
    @params['name'] = "Trinity_Assembly"
    @params['paired'] = ['false', 'true']
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['prinseqOpt'] = '-trim_qual_left 20 -trim_qual_right 20 -min_qual_mean 25 -min_len 50'
    @params['flexbarOpt'] = '--adapter-min-overlap 15 -at 1 --min-read-length 50'
    @params['trinityOpt'] = '--min_kmer_cov 2'
    @params['specialOptions'] = ''
    @params['mail'] = ""
 end
  def next_dataset
    {'Name'=>@params['name'],
     'Fasta [File]'=>File.join(@result_dir, "#{@params['name']}.fasta")
    }
  end
  def commands
    run_RApp
  end
end

if __FILE__ == $0
  usecase = EdgeRApp.new

  usecase.project = "p1001"
  usecase.user = 'masamasa'

  # set user parameter
  # for GUI sushi
  #usecase.params['process_mode'].value = 'SAMPLE'
  usecase.params['build'] = 'mm10'
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

