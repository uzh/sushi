#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class BamBoozleApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'BamBoozle'
    @params['process_mode'] = 'SAMPLE'
    @analysis_category = 'Prep'
    @description =<<-EOS
BAMboozle: Versatile removal of human sequence variation data for open data sharing  ToolLink: https://github.com/sandberg-lab/dataprivacy<br/>
    Reference:
    Ziegenhain, C., Sandberg, R. BAMboozle removes genetic variation from human sequence data for open data sharing. Nat Commun 12, 6216 (2021). https://doi.org/10.1038/s41467-021-26152-8 https://www.nature.com/articles/s41467-021-26152-8
    EOS
    @required_columns = ['Name', 'BAM', 'refBuild', 'Species', 'paired']
    @required_params = ['name', 'refBuild', 'paired']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '100'
    @params['name'] = 'BamBoozle'
    @params['paired'] = false
    @params['refBuild'] = ref_selector
    @params['cmdOptions'] = '--keepunmapped --strict'
    @params['cmdOptions', 'description'] = 'specify other commandline options for BAMBoozle'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Dev/R", "Dev/Python", "Tools/samtools", "Tools/BEDTools"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'Read1 [File]'=>File.join(@result_dir, "#{@dataset['Name']}_cleaned_R1.fastq.gz"),
     'Species'=>@dataset['Species'],
     'Read Count'=>@dataset['Read Count']
    }.merge(extract_columns(@inherit_tags))
    if @params['paired']
      dataset['Read2 [File]'] = File.join(@result_dir, "#{@dataset['Name']}_cleaned_R2.fastq.gz")
    end
    dataset
  end
  def set_default_parameters
     @params['refBuild'] = @dataset[0]['refBuild']
     @params['paired'] = @dataset[0]['paired']
  end
  def commands
    run_RApp("EzAppBamBoozle")
  end
end

if __FILE__ == $0

end
