#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-093844'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class CanuApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Canu'
    @analysis_category = 'Assemble'
    @description =<<-EOS
Canu long read genome assembler
<a href='http://canu.readthedocs.io/en/latest/quick-start.html'>http://canu.readthedocs.io/en/latest/quick-start.html</a>
EOS

    @required_columns = ['Name','Read1']
    @required_params = ['canuReadOpt', 'canuGenomeSize']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '50'
    @params['scratch'] = '400'
    @params['inputType'] = ["pacbioSmrtCell","fastqFile"]
    @params['inputType', 'description'] = 'direct output of a pacbio run or pre-processed (e.g., demultiplexed) fastq file?'
    @params['canuReadOpt'] = '-pacbio-raw'
    @params['canuReadOpt', 'description'] = 'input read types: -pacbio-raw, -pacbio-corrected, -nanopore-raw, -nanopore-corrected. Default is pacbio raw data'
    @params['canuGenomeSize'] = '5000'
    @params['canuGenomeSize', 'description'] = 'estimated genome size in Kbp'
    @params['cmdOptions'] = 'useGrid=false'
    @params['cmdOptions', 'description'] = 'specify other commandline options for Canu; do not specify any option that is already covered by the dedicated input fields'
    @params['mail'] = ""
    @modules = ["Assembly/Canu", "Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'Reads'=>@dataset['Reads'],
     'Draft [File]'=>File.join(@result_dir, "#{@dataset['Name']}.contigs.fasta"),
     'CanuOut [File]'=>File.join(@result_dir, "#{@dataset['Name']}"),
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppCanu")
  end
end

if __FILE__ == $0
  run CanuApp
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
