#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class DiffMethylationApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'DiffMethylation'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'DifferentialMethylation'
    @description =<<-EOS
    Finding differentially methylated regions via <a href='https://www.bioconductor.org/packages/release/bioc/html/dmrseq.html'/>Dmrseq</a> <br/>
    and <a href='https://www.bioconductor.org/packages/release/bioc/html/DSS.html'/>DSS</a> <br/>
    It requires the cov files from Bismark. <br/>
EOS
    @required_columns = ['Name','COV']
    @required_params = ['grouping', 'sampleGroup', 'refGroup', 'refBuild']
    @params['cores'] = '4'
    @params['ram'] = '80'
    @params['scratch'] = '100'
    @params['paired'] = true
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['grouping'] = ''
    @params['sampleGroup'] = ''
    @params['sampleGroup', 'description'] = 'sampleGroup should be different from refGroup'
    @params['refGroup'] = ''
    @params['refGroup', 'description'] = 'refGroup should be different from sampleGroup'
    @params['qVal'] = '0.1'
    @params['qVal', 'description'] = 'FDR cutoff for regions'
    @params['qVal_perSite'] = '0.001'
    @params['qVal_perSite', 'description'] = 'FDR cutoff for individual CpGs'
    @params['minDelta'] = '0.1'
    @params['minDelta', 'description'] = 'minimal difference to be reported'
    @params['cmdOptions'] = ''
    @params['mail'] = ""
    @modules = ["Dev/R"]
    @inherit_columns = ["Order Id"]
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
  end
  
  def next_dataset
    @comparison = "#{@params['sampleGroup']}--over--#{@params['refGroup']}"
    @params['comparison'] = @comparison
    @params['name'] = @comparison
    report_file = File.join(@result_dir, "#{@params['comparison']}")
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@params['name'],
      'Report [Link]'=>report_link,
      'ResultFolder [File]'=>report_file
    }.merge(extract_columns(colnames: @inherit_columns))
  end
  def commands
    run_RApp("EzAppDiffMethylation")
  end
end

if __FILE__ == $0
  
end
