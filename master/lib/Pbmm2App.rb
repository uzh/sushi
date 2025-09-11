#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class Pbmm2App < SushiFabric::SushiApp
  def initialize
    super
    @name = 'pbmm2'
    @analysis_category = 'Map'
    @description =<<-EOS
A minimap2 frontend for PacBio native data formats<br/>
<a href='https://github.com/PacificBiosciences/pbmm2'>Github website</a>
EOS

    @required_columns = ['Name','Read1']
    @required_params = ['refBuild','ReadOpt']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '100'
    @params['refBuild'] = ref_selector
    @params['refBuild', 'description'] = 'the genome refBuild and annotation to use as reference.'
    @params['ReadOpt'] = 'HIFI'
    @params['ReadOpt', 'description'] = 'input read types: SUBREAD, CCS, HIFI, ISOSEQ, UNROLLED. Default is HIFI'
    @params['cmdOptions'] = ""
    @params['mail'] = ""
    @modules = ["Tools/samtools"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def next_dataset
    {
     'Name'=>@dataset['Name'],
     'BAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam"),
     'BAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam.bai"),
     'IGV [File,Link]'=>File.join(@result_dir, "#{@dataset['Name']}-igv.html"),
     'refBuild'=>@params['refBuild'],
     'Pbmm2Log [File]'=>File.join(@result_dir, "#{@dataset['Name']}_pbmm2.log")
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppPbmm2")
  end
end

if __FILE__ == $0

end

