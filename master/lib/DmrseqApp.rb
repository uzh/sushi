#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class DmrseqApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'DmrseqApp'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'DifferentialMethylation'
    @description =<<-EOS
    Finding differentially methylated regions via <a href='https://www.bioconductor.org/packages/release/bioc/html/dmrseq.html'/>Dmrseq</a> <br/>
    It requires cov file from Bismark. <br/>
EOS
    @required_columns = ['Name','COV']
    @required_params = ['grouping', 'sampleGroup', 'refGroup']
    @params['cores'] = '4'
    @params['ram'] = '80'
    @params['scratch'] = '100'
    @params['paired'] = true
    @params['grouping'] = ''
    @params['sampleGroup'] = ''
    @params['sampleGroup', 'description'] = 'sampleGroup should be different from refGroup'
    @params['refGroup'] = ''
    @params['refGroup', 'description'] = 'refGroup should be different from sampleGroup'
    @params['qval'] = '0.05'
    @params['qval', 'description'] = 'FDR cutoff'
    @params['cmdOptions'] = ''
    @params['mail'] = ""
    @modules = ["Dev/R"]
    @inherit_columns = ["Order Id"]
  end
  def next_dataset
    @comparison = "#{@params['sampleGroup']}--over--#{@params['refGroup']}"
    @params['comparison'] = @comparison
    @params['name'] = @comparison
    {'Name'=>@params['name'],
      'CovObject [File]'=>'bismarkBSseq_filtered.rds',
      'RegionFile [File]'=>'dmrseq_results.rds',
      'BlockFile [File]'=>'large_blocks.rds'
    }.merge(extract_columns(colnames: @inherit_columns))
  end
  def commands
    run_RApp("EzAppDmrseq")
  end
end

if __FILE__ == $0
  
end
