#!/usr/bin/env ruby
# encoding: utf-8
Version = '20181127-155223'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class Samsa2App < SushiFabric::SushiApp
def initialize
super
@name = 'Samsa2'
@analysis_category = 'Metagenomics'
@description =<<-EOS
Metatranscriptomics pipeline with Samsa2.
<a href='https://github.com/transcript/samsa2'>Samsa2 on Github.</a>
EOS
@required_columns = ['Name', 'Read1']
@required_params = ['useSubsystemDB']
@params['cores'] = '8'
@params['cores', "context"] = "slurm"
@params['ram'] = '30'
@params['ram', "context"] = "slurm"
@params['scratch'] = '50'
@params['scratch', "context"] = "slurm"
@params['useSubsystemDB'] = false
@params['useSubsystemDB', 'description'] = 'Should the metatranscriptome be annotated also against Subsystem? Not supported at the moment.'
@params['useSubsystemDB', "context"] = "Samsa2"
@params['mail'] = ""
@inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
@modules = ["Dev/R", "Dev/jdk"]
end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  end
  def set_default_parameters
     @params['paired'] = dataset_has_column?('Read2')
  end
def next_dataset
   dataset = {'Name'=>@dataset['Name'],
     'annotationFileRefSeq [File]' => File.join(@result_dir, "#{@dataset['Name']}.RefSeq.annotated.txt"),
     'annotationORGFileRefSeq [File]' => File.join(@result_dir, "#{@dataset['Name']}.RefSeq.annotated_org.txt"),
     'annotationFUNCFileRefSeq [File]' => File.join(@result_dir, "#{@dataset['Name']}.RefSeq.annotated_func.txt")
}.merge(extract_columns(@inherit_tags))
  if @params['useSubsystemDB'] 
      dataset['annotationFileSubsystems [File]'] = File.join(@result_dir, "#{@dataset['Name']}.Subsys.annotated.txt")
  end
  dataset
end
def commands
run_RApp("EzAppSamsa2")
end
end

if __FILE__ == $0
end
