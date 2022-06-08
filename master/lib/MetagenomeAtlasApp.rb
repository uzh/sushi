#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class MetagenomeAtlasApp < SushiFabric::SushiApp
   def initialize
    super
    @name = 'MetagenomeAtlasApp'
    @analysis_category = 'Metagenomics'
    @description =<<-EOS
    Data preprocssing with Metagenome Atlas. The input is short reads, requires both R1 and R2 reads and assembles genomes found in the sample.
    <a href='https://github.com/metagenome-atlas/atlas'>Metagenome Atlas main repository.</a>
EOS
    @params['process_mode'] = 'DATASET'
    @required_columns = ['Name', 'Read1', 'Read2']
    @required_params = ['paired']    
    @params['cores'] = '8'
    @params['paired'] = 'True'
    @params['paired', 'description'] = 'whether the reads are paired end; if false then only Read1 is considered even if Read2 is available.'
    @params['ram'] = '60'
    @params['scratch'] = '60'
    @params['mail'] = ""
    @params['name'] = "MetagenomeAtlas"
    @inherit_tags = ['B-Fabric']
    @modules = ['Dev/R', 'QC/mash/2.3' , 'QC/fastANI/1.32']
  end

def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@params['name'],
     'Static Report [Link]'=>report_link,
     'Report [File]'=>report_file,
    }    
end
  def commands
     run_RApp("EzAppMetagenomeAtlas", conda_env: "metagenome-atlas")
  end
end

if __FILE__ == $0

end


