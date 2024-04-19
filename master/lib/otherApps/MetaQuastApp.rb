#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-095008'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class MetaQuastApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'MetaQuast'
    @analysis_category = 'QC'
    @description =<<-EOS
MetaQUAST evaluates and compares metagenome assemblies based on alignments to close references
<a href='https://quast.sourceforge.net/metaquast.html'>https://quast.sourceforge.net/metaquast.html</a>
EOS
    @required_columns = ['Name','Draft']
    @required_params = ['cores', 'ram', 'scratch']
    # optional params
    @params['cores'] = '4'
    @params['ram'] = '30'
    @params['scratch'] = '50'
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'specify other commandline options for QUAST; do not specify any option that is already covered by the dedicated input fields'
    @params['mail'] = ""
    @modules = ["Aligner/ncbi-blast/2.12.0+","QC/QUAST/5.2.0", "Dev/R", "Dev/Python/2.7.16"]
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
    run_RApp("EzAppMetaQuast")
  end
end

if __FILE__ == $0
  run MetaQuastApp
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
