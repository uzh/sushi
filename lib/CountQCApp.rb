#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class CountQCApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'CountQC'
    @analysis_category = 'QC'
    @params['process_mode'] = 'DATASET'
    @required_columns = ['Name','Count', 'Species', 'build', 'featureLevel', 'featureFile']
    @required_params = []
    # optional params
    @params['cores'] = '1'
    @params['ram'] = '2'
    @params['scratch'] = '10'
    @params['name'] = 'Count_QC'
    @params['build'] = ref_selector
    @params['featureFile'] = 'genes.gtf'
    @params['featureLevel'] = ['gene', 'isoform']
    @params['normMethod'] = 'logMean'
    @params['expressionName'] = ''
    @params['runGO'] = false
    @params['specialOptions'] = ''
    @params['mail'] = ""
  end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@params['name'],
     'Species'=>@dataset['Species'],
     'build'=>@params['build'],
     'Report [File]'=>report_file,
     'Html [Link]'=>report_link,
    }
  end
  def set_default_parameters
    @params['build'] = @dataset[0]['build']
    if dataset_has_column?('featureFile')
      @params['featureFile'] = @dataset[0]['featureFile']
    end
  end
  def commands
    command = "/usr/local/ngseq/bin/R --vanilla --slave << EOT\n"
    command << "R_SCRIPT_DIR <<- '#{GlobalVariables::R_SCRIPT_DIR}'\n"
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
    command << "tryCatch({countQCApp(input=inputDatasetFile, output=output, config=config)}, error=function(e){my.mail(to=config[['mail']], text=e); stop(e)})\n"
    command << "EOT"
    command
  end
end

if __FILE__ == $0

end

