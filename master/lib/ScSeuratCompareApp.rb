#!/usr/bin/env ruby
# encoding: utf-8
Version = '20210702-145619'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class ScSeuratCompareApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'ScSeuratCompare'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
    Empirical analysis of digital gene expression data in R<br/>
<a href='https://bioconductor.org/packages/release/bioc/html/DifferentialState.html'>DifferentialState</a><br/>
    EOS
    @required_columns = ['Name', 'Report']
    @required_params = ['CellIdentity', 'grouping', 'sampleGroup', 'refGroup', 'pseudoBulkMode']
    # optional params
    @params['cores'] = '4'
    @params['ram'] = '30'
    @params['scratch'] = '50'
    @params['DE.method'] = ['wilcox', 'LR']
    @params['DE.method', 'description'] = "Method to be used when calculating differentially expressed genes between conditions."
    @params['DE.regress'] = ['Batch', 'CellCycle']
    @params['DE.regress','multi_selection'] = true
    @params['DE.regress', 'description'] = "Variables to regress when calculating differentially expressed genes. Only used with the LR method."
    @params['CellIdentity'] = 'ident'
    @params['CellIdentity', 'description'] = "The Seurat metadata column which contains the cell clusters or types (usually or 'ident' 'cellTypes' if the clusters have been labeled)"
    @params['grouping'] = 'Condition'
    @params['grouping', 'description'] = "The Seurat metadata column which contains the sample grouping information"
    @params['pseudoBulkMode'] = false
    @params['pseudoBulkMode', 'description'] = "Weather to aggregate the counts to the pseudo-bulk level prior to performing the DE experiments. Setting this to true also requires setting 'replicateGrouping'"
    @params['replicateGrouping'] = ""
    @params['replicateGrouping', 'description'] = "(pseudo-bulk mode only) The column in the Seurat metadata containing the replicate group information."
    @params['sampleGroup'] = ''
    @params['sampleGroup', 'description'] = 'sampleGroup should be different from refGroup'
    @params['refGroup'] = ''
    @params['refGroup', 'description'] = 'refGroup should be different from sampleGroup'
    @params['mail'] = ""
    @params['Rversion'] = ["Dev/R/4.4.0", "Dev/R/4.3.0", "Dev/R/4.2.2", "Dev/R/4.2.0", "Dev/R/4.1.2", "Dev/R/4.1.0"]
    @inherit_columns = ["Order Id"]
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
    }.merge(extract_columns(colnames: @inherit_columns))
  end
  def commands
    command = "module load #{@params["Rversion"]}\n"
    command << run_RApp("EzAppScSeuratCompare")
  end
end

if __FILE__ == $0

end
