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
    @analysis_category = 'Metagenomics'
    @description =<<-EOS
MetaQUAST (Quality Assessment Tool for Metagenome Assemblies)
<a href='http://quast.sourceforge.net/metaquast'>http://quast.sourceforge.net/metaquast</a>
EOS
    @required_columns = ['Name','contigFile']
    @required_params = ['cores', 'ram', 'scratch']
    @params['cores'] = '4'
    @params['ram'] = '30'
    @params['scratch'] = '50'
    @params['fileWithListOfRefs'] = ''
    @params['fileWithListOfRefs', 'description'] = 'full path to a txt file listing the reference genomes.'
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'specify other commandline options for Metaquast; do not specify any option that is already covered by the dedicated input fields'
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
    run_RApp("EzAppMetaQuast")
  end
end

if __FILE__ == $0
end
