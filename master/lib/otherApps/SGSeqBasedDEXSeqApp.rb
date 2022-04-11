#!/usr/bin/env ruby
# encoding: utf-8
Version = '20180412-100155'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class SGSeqBasedDEXSeqApp < SushiFabric::SushiApp
def initialize
super
@name = 'SGSeqBasedDEXSeq'
@analysis_category = 'Alternative_Splicing'
@description =<<-EOS
Inference of differential exon usage in RNA-Seq<br/>
    <a href='https://bioconductor.org/packages/release/bioc/html/DEXSeq.html'>DEXSeq</a><br/>
    EOS
@params['process_mode'] = 'DATASET'
@required_columns = ['Name','SgSeqRDataFile','SgSeqCountFile']
@required_params = ['grouping', 'sampleGroup', 'refGroup','refBuild']
@params['cores'] = '4'
@params['ram'] = '30'
@params['scratch'] = '100'
@params['refBuild'] = ref_selector
@params['refFeatureFile'] = 'genes.gtf'
@params['paired'] = false
@params['strandMode'] = ['both', 'sense', 'antisense']
@params['grouping'] = ''
@params['sampleGroup'] = ''
@params['sampleGroup', 'description'] = 'sampleGroup should be different from refGroup'
@params['refGroup'] = ''
@params['refGroup', 'description'] = 'refGroup should be different from sampleGroup'
@params['mail'] = ""
@modules = ["Dev/R"]
end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('refFeatureFile')
      @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    end
    if dataset_has_column?('paired')
      @params['paired'] = @dataset[0]['paired']
    end
    if dataset_has_column?('strandMode')
      @params['strandMode'] = @dataset[0]['strandMode']
    end
  end
  
def next_dataset
@comparison = "#{@params['sampleGroup']}--over--#{@params['refGroup']}"
@params['comparison'] = @comparison
@params['name'] = @comparison
report_file = File.join(@result_dir, "#{@params['comparison']}")
report_link = File.join(report_file, '00index.html')
{'Name'=>@comparison,
    'DEXSeqResults [File]'=>File.join(@result_dir, "#{@comparison}.DEXSeqResults.txt"),
}
end
def commands
run_RApp("EzAppSGseqBasedDEXSeq")
end
end

if __FILE__ == $0

end
