#!/usr/bin/env ruby
# encoding: utf-8
Version = '20210702-145619'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class ScCompareSCEApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'ScCompareSCE'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Differential_Expression'
    @description =<<-EOS
    Empirical analysis of digital gene expression data in R<br/>
<a href='https://bioconductor.org/packages/release/bioc/html/DifferentialState.html'>DifferentialState</a><br/>
    EOS
    @required_columns = ['Name', 'Condition', 'Batch']
    @required_params = ['grouping', 'sampleGroup', 'refGroup']
    # optional params
    @params['cores'] = '4'
    @params['ram'] = '30'
    @params['scratch'] = '50'
    @params['grouping'] = '' ## Note: this will be a selector defined by Factor tagged column
    @params['sampleGroup'] = '' ## Note: this will be a selector defined by Factor tagged column
    @params['sampleGroup', 'description'] = 'sampleGroup should be different from refGroup'
    @params['refGroup'] = '' ## Note: this will be a selector defined by Factor tagged column
    @params['refGroup', 'description'] = 'refGroup should be different from sampleGroup'
    @params['mail'] = ""
    @params['Rversion'] = ["Dev/R/4.1.2", "Dev/R/4.1.0", "Dev/R/4.0.4", "Dev/R/4.0.3"]
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
    command << run_RApp("EzAppScCompareSCE")
  end
end

if __FILE__ == $0

end
