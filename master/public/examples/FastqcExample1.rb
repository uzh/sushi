#!/usr/bin/env ruby
# encoding: utf-8
Version = '20160425-055647'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class FastqcExample1 < SushiFabric::SushiApp
  def initialize
    super
    @name = 'FastqcExample1'
    @analysis_category = 'Map'
    @description =<<-EOS

    EOS
    @required_columns = ['Name','Read1']
    @required_params = ['paired']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '15'
    @params['scratch'] = '100'
    @params['library_type'] = ['fr-unstranded', 'fr-firststrand', 'fr-secondstrand']
    @params['paired'] = false
    @params['paired', 'description'] = 'whether the reads are paired-end or single-end'
  end
  def set_default_parameters
    @params['paired'] = dataset_has_column?('Read2')
  end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  end
  def next_dataset
    {'Name'=>@dataset['Name'], 
     'Report [Link]'=>File.join(@result_dir, "#{@dataset['Name']}/fastqc_report.html"), 
     'ReportDir [File]'=>File.join(@result_dir, "#{@dataset['Name']}")
    }
  end
  def commands
    cmd = "fastqc --extract -o . -t #{@params['cores']} #{@dataset['Read1']}"
    cmd << "\n"
    fastqcDir = File.basename(@dataset['Read1'].to_s).gsub('.fastq.gz','_fastqc')
    cmd << "mv #{fastqcDir} #{@dataset['Name']}"
    cmd << "\n"
  end
end
