#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class NcPROApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'ncPRO_Report'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Count'
    @description =<<-EOS
Annotation and Profiling of ncRNAs in smallRNA-seq<br/>
Uses ncPRO-seq for a complete analysis of small-RNA-seq. ncPRO-seq considers performs quality assessment and quantiation.
It considers various species of small RNA.<br />
<a href='https://ncpro.curie.fr/index.html'>https://ncpro.curie.fr/index.html/</a>
ncPRO can only process single-end stranded RNA.
EOS

    @required_columns = ['Name','Read1', 'Adapter1', 'Species']
    @required_params = ['name', 'cores', 'ram', 'scratch']

    @params['cores'] = '8'
    @params['ram'] = '16'
    @params['scratch'] = '100'
    @params['build'] = ref_selector
    #@params['paired'] = false
    #@params['paired', 'description'] = 'whether the reads are paired end; must be false since ncPRO-seq does not support paired-end'
    @params['name'] = 'ncPRO_Result'
    @params['mail'] = ""
  end
  def preprocess
    if @params['paired']
      @required_columns<<  'Read2'
    end
  end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, 'ncpro/report.html')
    {'Name'=>@params['name'],
     'Species'=>@dataset['Species'],
     'build'=>@params['build'],
     'Report [File]'=>report_file,
     'Html [Link]'=>report_link,
     'TrimCounts [Link]'=>File.join(report_file, 'trimCounts-barplot.png')
    }
  end
  def commands
    command = "/usr/local/ngseq/bin/R --vanilla --slave<<  EOT\n"
    command << "R_SCRIPT_DIR <<- '#{R_SCRIPT_DIR}'\n"
    command<<  "source(file.path(R_SCRIPT_DIR, 'init.R'))\n"
    command << "config = list()\n"
    config = @params
    config.keys.each do |key|
      command << "config[['#{key}']] = '#{config[key]}'\n" 
    end
    command << "config[['dataRoot']] = '#{@gstore_dir}'\n"
    command << "output = list()\n"
    output = next_dataset
    output.keys.each do |key|
      command << "output[['#{key}']] = '#{output[key]}'\n" 
    end
    command<<  "inputDatasetFile = '#{@input_dataset_tsv_path}'\n"
    command<<  "ncproApp(input=inputDatasetFile, output=output, config=config)\n"
    command<<  "EOT\n"
    command
  end
end

if __FILE__ == $0
  run NcPROApp

end




