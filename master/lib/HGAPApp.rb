#!/usr/bin/env ruby
# encoding: utf-8
Version = '20170311'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class HGAPApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'HGAP'
    @analysis_category = 'Assemble'
    @description =<<-EOS
Canu long read genome assembler
<a href='http://canu.readthedocs.io/en/latest/quick-start.html'>http://canu.readthedocs.io/en/latest/quick-start.html</a>
EOS
    
    @required_columns = ['Name','Reads']
    @required_params = ['genomeSize', 'xCoverage']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '50'
    @params['scratch'] = '400'
    @params['genomeSize'] = '5000000'
    @params['genomeSize', 'description'] = 'The approximate genome size, in base pairs.' 
    @params['xCoverage'] = '25'
    @params['xCoverage', 'description'] = 'Fold coverage to target for when picking the minimum fragment length for assembly; typically 15 to 25.'
    @params['mail'] = ""
  end
  def next_dataset
    {'Name'=>@dataset['Name'], 
     'Reads'=>@dataset['Reads'], 
     'Draft [File]'=>File.join(@result_dir, "#{@dataset['Name']}", "data", "polished_assembly.fasta.gz"),	
     'HGAPOut [Html]'=>File.join(@result_dir, "#{@dataset['Name']}"),
    }.merge(extract_column("Factor")).merge(extract_column("B-Fabric"))
  end
  def commands
    run_RApp("EzAppHGAP")
  end
end

if __FILE__ == $0
  run HGAPApp
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

