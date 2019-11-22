#!/usr/bin/env ruby
# encoding: utf-8
Version = '20180209-163825'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class SGSeqApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'SGSeq'
    @analysis_category = 'Alternative_Splicing'
 @description =<<-EOS
    Data preparation for alternative splicing analysis<br/>
<a href='https://bioconductor.org/packages/devel/bioc/vignettes/SGSeq/inst/doc/SGSeq.html'>manual</a><br/>
EOS
    @required_columns = ['Name','BAM','BAI']
    @required_params = ['refBuild','refFeatureFile']
    @params['cores'] = '8'
    @params['ram'] = '20'
    @params['scratch'] = '10'
    @params['refBuild'] = ref_selector
    @params['paired'] = false
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['refFeatureFile'] = 'genes.gtf'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end

  def next_dataset
    {'Name'=>@dataset['Name'],
     'SgSeqRDataFile [File]'=>File.join(@result_dir, "#{@dataset['Name']}.Rdata"),
     'SgSeqVarFreqFile [File]'=>File.join(@result_dir, "#{@dataset['Name']}.varFreq.txt"),
     'SgSeqCountFile [File]'=>File.join(@result_dir, "#{@dataset['Name']}.inputForDEXSeq.txt"),
     'refBuild'=>@params['refBuild'],
     'refFeatureFile'=>@params['refFeatureFile']
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppSGSeq")
  end
end

if __FILE__ == $0

end
