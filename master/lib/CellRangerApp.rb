#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class CellRangerApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'CellRangerCount'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
This wrapper runs <a href='https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count',>cellranger count</a> in Single-library analysis mode.
    EOS
    @required_columns = ['Name','RawDataDir','Species']
    @required_params = ['name']
    @params['cores'] = '8'
    @params['ram'] = '32'
    @params['scratch'] = '100'
    @params['name'] = 'CellRangerCount'
    @params['reference'] = {'select'=>''}
    Dir["/srv/GT/databases/10X_References/*"].sort.select{|10X| File.directory?(10X)}.each do |dir|
      @params['reference'][File.basename(dir)] = File.basename(dir)
    end
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Dev/R"]
    @inherit_tags = ["B-Fabric"]
  end
  def set_default_parameters
  end
  def next_dataset
    report_dir = File.join(@result_dir,"#{@dataset['Name']}")
    dataset = {
      'Name'=>@dataset['Name'],
      'ResultDir [File]'=>report_dir,
      'Report [Link]'=>File.join(report_dir, 'outs/web_summary.html'),
      'BAM [File]'=>File.join(report_dir, 'outs/possorted_genome_bam.bam'),
      'BAI [File]'=>File.join(report_dir, 'outs/possorted_genome_bam.bam.bai'),
      'Species'=>@dataset['Species'],
      'reference'=>@params['reference']
    }.merge(extract_columns(@inherit_tags))
    dataset
  end
  def commands
    run_RApp("EzAppCellRanger")
  end
end

if __FILE__ == $0

end
