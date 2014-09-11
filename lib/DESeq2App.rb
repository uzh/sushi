#!/usr/bin/env ruby
# encoding: utf-8
Version = '20131128-084558'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class DESeq2App < SushiFabric::SushiApp
  def initialize
    super
    @name = 'DESeq2'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Differential_Expression'
    @required_columns = ['Name','Count', 'Species', 'build', 'featureLevel', 'featureFile']
    @required_params = ['grouping', 'sampleGroup', 'refGroup']
    # optional params
    @params['cores'] = '1'
    @params['ram'] = '2'
    @params['scratch'] = '10'
    @params['build'] = ref_selector
    @params['featureFile'] = 'genes.gtf'
    @params['featureLevel'] = ['gene', 'isoform']
    @params['grouping'] = '' 
    @params['sampleGroup'] = '' 
    @params['refGroup'] = ''
    #@params['normMethod'] = 'logMean'
    @params['runGO'] = ['false', 'true']
    @params['specialOptions'] = ''
    @params['mail'] = ""
  end
  def next_dataset
    @comparison = "#{@params['sampleGroup']}--over--#{@params['refGroup']}"
    @params['comparison'] = @comparison
    report_file = File.join(@result_dir, @comparison)
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@comparison,
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
    command << "source('/usr/local/ngseq/opt/sushi_scripts/init.R')\n"
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
    command << "DESeq2App(input=inputDatasetFile, output=output, config=config)\n"
    command << "EOT"
    command
  end
end

if __FILE__ == $0

end

