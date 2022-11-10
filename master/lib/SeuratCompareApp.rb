#!/usr/bin/env ruby
# encoding: utf-8
Version = '20210702-145619'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class SeuratCompareApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'SeuratCompare'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
    Empirical analysis of digital gene expression data in R<br/>
<a href='https://bioconductor.org/packages/release/bioc/html/DifferentialState.html'>DifferentialState</a><br/>
    EOS
    @required_columns = ['Name', 'Report']
    @required_params = ['grouping', 'sampleGroup', 'refGroup']
    # optional params
    @params['cores'] = '4'
    @params['ram'] = '30'
    @params['scratch'] = '50'
    @params['DE.method'] = ['wilcox', 'LR']
    @params['DE.method', 'description'] = "Method to be used when calculating differentially expressed genes between conditions."
    @params['DE.regress'] = ['Batch', 'CellCycle']
    @params['DE.regress','multi_selection'] = true
    @params['DE.regress', 'description'] = "Variables to regress when calculating differentially expressed genes. Only used with the LR method."
    @params['grouping'] = '' ## Note: this will be a selector defined by Factor tagged column
    @params['sampleGroup'] = '' ## Note: this will be a selector defined by Factor tagged column
    @params['sampleGroup', 'description'] = 'sampleGroup should be different from refGroup'
    @params['refGroup'] = '' ## Note: this will be a selector defined by Factor tagged column
    @params['refGroup', 'description'] = 'refGroup should be different from sampleGroup'
    @params['mail'] = ""
    @params['Rversion'] = ["Dev/R/4.2.0", "Dev/R/4.1.2", "Dev/R/4.1.0"]
  end
   def preprocess
    @random_string = (1..12).map{[*('a'..'z')].sample}.join
  end
  def next_dataset
    @comparison = "#{@params['sampleGroup']}--over--#{@params['refGroup']}"
    @params['comparison'] = @comparison
    @params['name'] = @comparison
    report_file = File.join(@result_dir, "#{@params['comparison']}")
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@comparison,
     'Static Report [Link]'=>report_link,
     'Report [File]'=>report_file,
    }
  end
  def commands
    command = "module load #{@params["Rversion"]}\n"
    command << run_RApp("EzAppSeuratCompare")
  end
end

if __FILE__ == $0

end
