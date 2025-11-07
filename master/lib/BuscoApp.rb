#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-095008'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class BuscoApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Busco'
    @analysis_category = 'QC'
    @description =<<-EOS
BUSCO: from QC to gene prediction and phylogenomics
<a href='https://busco.ezlab.org/'>https://busco.ezlab.org/</a>
EOS
    @required_columns = ['Name','Draft']
    @required_params = ['cores', 'ram', 'scratch', 'mode', 'lineage']
    # optional params
    @params['cores'] = '4'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '30'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '50'
    @params['scratch', "context"] = "slurm"
    @params['mode'] = ['geno', 'tran',  'prot']
    @params['mode', 'description'] = 'analysis mode: genome, transcriptome, proteins'
    @params['mode', "context"] = "Busco"
    @params['lineage'] = ''
    @params['lineage', 'description'] = 'the name of the BUSCO lineage dataset to be used, such as bacteria_odb10 or more specifically gammaproteobacteria_odb10. A complete list of available datasets can be found under https://fgcz-gstore.uzh.ch/reference/BUSCO/v5/busco_datasets.txt'
    @params['lineage', "context"] = "Busco"
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'specify other commandline options for prokka; do not specify any option that is already covered by the dedicated input fields'
    @params['cmdOptions', "context"] = "Busco"
    @params['mail'] = ""
    @modules = ["Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'BuscoPlot [Link]'=>File.join(@result_dir, "#{@dataset['Name']}", 'busco_figure.png'),
     'BuscoOut [File]'=>File.join(@result_dir, "#{@dataset['Name']}"),
     'Species'=>@dataset['Species'],
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppBusco",conda_env: "gi_busco5.8.0")
  end
end

if __FILE__ == $0
  run BuscoApp
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
