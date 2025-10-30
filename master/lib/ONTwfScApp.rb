#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class ONTwfScApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'ONTwfSc'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
A research pipeline designed to identify the cell barcode and UMI sequences present in nanopore sequencing reads generated from single-cell gene expression libraries<br/>
<a href='https://github.com/epi2me-labs/wf-single-cell'>ONT Workflow single-cell</a>
EOS

    @required_columns = ['Name','Read1']
    @required_params = ['refBuild']
    # optional params
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '30'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '100'
    @params['scratch', "context"] = "slurm"
    @params['refBuild'] = 'refdata-gex-GRCh38-2020-A'
    @params['refBuild', 'description'] = 'the genome refBuild and annotation to use as reference.'
    @params['refBuild', "context"] = "referfence genome assembly"
    @params['kitName'] = ['3prime', '5prime', 'multiome']
    @params['kitName', 'description'] = '10x kit name'
    @params['kitName', "context"] = "ONTwfSc"
    @params['kitVersion'] = ['v3', 'v2', 'v1']
    @params['kitVersion', 'description'] = '10x kit version'
    @params['kitVersion', "context"] = "ONTwfSc"
    @params['expCells'] = '500'
    @params['expCells', 'description'] = 'Number of expected cells.'
    @params['expCells', "context"] = "ONTwfSc"
    @params['cmdOptions'] = ""
    @params['cmdOptions', "context"] = "ONTwfSc"
    @params['mail'] = ""
    @modules = ["Dev/jdk", "Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def next_dataset
    {
    'Name'=>@dataset['Name'],
     'OutDir [File]'=>File.join(@result_dir, "#{@dataset['Name']}"),
     'OutReport [Link]'=>File.join(@result_dir, "#{@dataset['Name']}", "wf-single-cell-report.html")
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppONTwfSc")
  end
end

if __FILE__ == $0

end

