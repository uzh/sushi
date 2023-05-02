#!/usr/bin/env ruby
# encoding: utf-8
Version = '20180718-144150'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class CNVnatorApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'CNVnator'
    @analysis_category = 'Variants'
    @description =<<-EOS
 A tool for CNV discovery and genotyping from depth-of-coverage by mapped reads <br/>
<a href='https://github.com/abyzovlab/CNVnator'>https://github.com/abyzovlab/CNVnator</a>
    EOS
    @required_columns = ['Name','BAM','BAI', 'refBuild']
    @required_params = ['refBuild','paired']
    # optional params
    @params['cores'] = '1'
    @params['ram'] = '15'
    @params['scratch'] = '100'
    @params['refBuild'] = ref_selector
    @params['paired'] = false
    @params['refFeatureFile'] = 'genes.gtf'
    @params['binSize'] = 1000
    @params['maxEVal'] = 0.01
    @params['maxQ0'] = 0.1
    @params['cmdOptions'] = ''
    @params['specialOptions'] = ''
    @params['mail'] = ''
    @modules = ["Dev/R"]

    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('refFeatureFile')
      @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    end
    if dataset_has_column?('paired')
      @params['paired'] = @dataset[0]['paired']
    end
 end

  def next_dataset
    filteredCNV_link = File.join(@result_dir, "#{@dataset['Name']}_CNV.txt")
    rawCNV_link = File.join(@result_dir, "#{@dataset['Name']}.txt")

    {'Name'=>@dataset['Name'],
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'refFeatureFile'=>@params['refFeatureFile'],
     'paired'=>@params['paired'],
     'CalledCNVs [File]'=>rawCNV_link,
     'FilteredCNVs [File]'=>filteredCNV_link,
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    command =<<-EOS 
    #source /usr/local/ngseq/src/root_v6.26.10/bin/thisroot.sh
    #export YEPPPLIBDIR=/usr/local/ngseq/src/yeppp-1.0.0/binaries/linux/x86_64
    #export LD_LIBRARY_PATH=$YEPPPLIBDIR:$LD_LIBRARY_PATH 
    EOS
    command + run_RApp("EzAppCNVnator")
  end
end

if __FILE__ == $0

end
