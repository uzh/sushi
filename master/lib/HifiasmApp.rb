#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-093844'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class HifiasmApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Hifiasm'
    @analysis_category = 'Assemble'
    @description =<<-EOS
Hifiasm long read genome assembler
<a href='https://hifiasm.readthedocs.io/en/latest/index.html'>https://hifiasm.readthedocs.io/en/latest/index.html</a>
EOS

    @required_columns = ['Name','Read1']
    @required_params = ['inputType']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '400'
    @params['inputType'] = ["HiFi", "ONT"]
    @params['inputType', 'description'] = 'PacBio HiFi reads or Oxford Nanopore reads'
    @params['ploidy'] = '2'
    @params['ploidy', 'description'] = 'number of haplotypes'
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'specify other commandline options for hifiasm; do not specify any option that is already covered by the dedicated input fields'
    @params['mail'] = ""
    @modules = ["Assembly/hifiasm", "Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'Reads'=>@dataset['Reads'],
     'Draft [File]'=>File.join(@result_dir, "#{@dataset['Name']}.p_ctg.fa"),
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppHifiasm")
  end
end

if __FILE__ == $0
  run HifiasmApp
  #usecase = Bowtie2App.new

  #usecase.project = "p1001"
  #usecase.user = 'masamasa'

  # set user parameter
  # for GUI sushi
  #usecase.params['process_mode'].value = 'SAMPLE'
  #usecase.params['refBuild'] = 'mm10'
  #usecase.params['paired'] = true
  #usecase.params['strandMode'] = 'both'
  #usecase.params['cores'] = 8
  #usecase.params['node'] = 'fgcz-c-048'

  # also possible to load a parameterset csv file
  # mainly for CUI sushi
  #usecase.parameterset_tsv_file = 'tophat_parameterset.tsv'
  #usecase.parameterset_tsv_file = 'test.tsv'

  # set input dataset
  # mainly for CUI sushi
  #usecase.dataset_tsv_file = 'tophat_dataset.tsv'

  # also possible to load a input dataset from Sushi DB
  #usecase.dataset_sushi_id = 3

  # run (submit to workflow_manager)
  #usecase.run
  #usecase.test_run

end
