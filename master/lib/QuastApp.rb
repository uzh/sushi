#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-095008'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class QuastApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Quast'
    @analysis_category = 'QC'
    @description =<<-EOS
QUAST (Quality Assessment Tool for Genome Assemblies)
<a href='http://quast.bioinf.spbau.ru/manual.html '>http://quast.bioinf.spbau.ru/manual.html</a>
EOS
    @required_columns = ['Name','Draft']
    @required_params = ['cores', 'ram', 'scratch']
    # optional params
    @params['cores'] = '4'
    @params['ram'] = '30'
    @params['scratch'] = '50'
    @params['refGenome'] = ''
    @params['refGenome', 'description'] = 'full path to a reference genome as a multi-fasta file'
    @params['refGene'] = ''
    @params['refGene', 'description'] = 'full path to a gene annotation file of the reference genome. Must be in gff or bed format'
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'specify other commandline options for QUAST; do not specify any option that is already covered by the dedicated input fields'
    @params['mail'] = ""
    @modules = ["QC/QUAST", "Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'QuastReport [Link]'=>File.join(@result_dir, "#{@dataset['Name']}", 'report.html'),
     'QuastOut [File]'=>File.join(@result_dir, "#{@dataset['Name']}"),
     'Species'=>@dataset['Species'],
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppQuast")
  end
end

if __FILE__ == $0
  run QuastApp
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
