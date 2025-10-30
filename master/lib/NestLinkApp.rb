#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class NestLinkApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'NestLinkApp'
    @analysis_category = 'Misc'
    @description =<<-EOS
NestLink - an R data package to guide through Engineered Peptide Barcodes for In-Depth Analyzes of Binding Protein Ensembles - https://bioconductor.org/packages/release/data/experiment/html/NestLink.html
    EOS
    @required_columns = ['Name','Read1', 'FlashLog']
    @required_params = ['NB_Linker1', 'NB_Linker2', 'ProteaseSite','FC_Linker','knownNBPath']
    @params['cores'] = '1'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '40'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '30'
    @params['scratch', "context"] = "slurm"
    @params['nReads'] = '0'
    @params['nReads', "context"] = "NestLink"
    @params['maxMismatch'] = '1'
    @params['maxMismatch', "context"] = "NestLink"
    @params['NB_Linker1'] = "GGCCggcggGGCC"
    @params['NB_Linker2'] = "GCAGGAGGA"
    @params['ProteaseSite'] = "TTAGTCCCAAGA"
    @params['FC_Linker'] = "GGCCaaggaggcCGG"
    @params['knownNBPath'] = '/home/huberlea/knownNB_new_corrected_LO.txt'
    @params['minRelBestHitFreq'] = '0'
    @params['minConsensusScore'] = '0.9'
    @params['minNanobodyLength'] = '321'
    @params['minFlycodeLength'] = '33'
    @params['FCminFreq'] = '4'
    @params['name'] = 'NestLink'
    @params['specialOptions'] = ''
    @params['mail'] = ''
    @modules = ["Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'NestLink Result [File]'=>File.join(@result_dir, "#{@dataset['Name']}")
    }.merge(extract_columns(@inherit_tags))
  end
  def set_default_parameters
  end
  def commands
    run_RApp("EzAppNestLink")
  end
end

if __FILE__ == $0

end
