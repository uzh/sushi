#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class ONTwfAAVqcApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'ONTwfAAVqc'
    @analysis_category = 'QC'
    @description =<<-EOS
QC of recombinant adeno-associated viral vector (rAAV) preparations.<br/>
<a href='https://nanoporetech.com/document/epi2me-workflows/wf-aav-qc'>wf-aav-qc</a>
EOS

    @required_columns = ['Name','Read1']
    @required_params = ['refHelper','refRepCap','refTrans','refBuild','itr1Start','itr1End','itr2Start','itr2End']
    # optional params
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '32'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '100'
    @params['scratch', "context"] = "slurm"
    @params['refHelper'] = '/srv/GT/databases/wf-aav-qc/helper.fasta'
    @params['refHelper', 'description'] = 'The helper plasmid FASTA file.'
    @params['refRepCap'] = '/srv/GT/databases/wf-aav-qc/repcap.fasta'
    @params['refRepCap', 'description'] = 'The rep/cap plasmid FASTA file.'
    @params['refTrans'] = '/srv/GT/databases/wf-aav-qc/transgene.fasta'
    @params['refTrans', 'description'] = 'The transgene plasmid FASTA file.'
    @params['refBuild'] = ref_selector
    @params['refBuild', 'description'] = 'Host cell reference genome.'
    @params['refBuild', "context"] = "reference genome assembly"
    @params['itr1Start'] = '11'
    @params['itr1Start', 'description'] = 'The start position of ITR1.'
    @params['itr1End'] = '156'
    @params['itr1End', 'description'] = 'The end position of ITR1.'
    @params['itr2Start'] = '2156'
    @params['itr2Start', 'description'] = 'The start position of ITR2.'
    @params['itr2End'] = '2286'
    @params['itr2End', 'description'] = 'The end position of ITR2.'
    @params['cmdOptions'] = ""
    @params['cmdOptions', "context"] = "ONTwfAAVqc"
    @params['mail'] = ""
    @modules = ["Dev/jdk", "Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def next_dataset
    {
    'Name'=>@dataset['Name'],
     'OutDir [File]'=>File.join(@result_dir, "#{@dataset['Name']}"),
     'OutReport [Link]'=>File.join(@result_dir, "#{@dataset['Name']}", "wf-aav-qc-report.html")
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppONTwfAAVqc")
  end
end

if __FILE__ == $0

end

