#!/usr/bin/env ruby
# encoding: utf-8
Version = '20140317-150018'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class PileupApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Pileup'
    @analysis_category = 'Variant_Analysis'
    @description =<<-EOS
Does simple variant calling using samtools mpileup and bcftools. Performs no variant annotation.
EOS
    
    @required_columns = ['Name','BAM','BAI', 'build']
    @required_params = []
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '16'
    @params['scratch'] = '100'
    @params['build'] = ref_selector
    @params['build', 'description'] = 'the genome build and annotation to use as reference'
    @params['paired'] = false
    @params['paired', 'description'] = 'whether the reads are paired end; if false then only Read1 is considered even if Read2 is available.'
    @params['mpileupOptions'] = ''
    @params['bcftoolsOptions'] = '-c'
    @params['rmdup'] = true
    @params['rmdup', 'description'] = 'remove duplicates?'
    @params['specialOptions'] = ''
    @params['specialOptions', 'description'] = 'special unsupported options that the R wrapper may support, format: <key>=<value>'
  end
 def set_default_parameters
     @params['build'] = @dataset[0]['build']
    if dataset_has_column?('featureFile')
      @params['featureFile'] = @dataset[0]['featureFile']
    end
    if dataset_has_column?('paired')
      @params['paired'] = @dataset[0]['paired']
    end                               
  end
  def next_dataset
    {'Name'=>@dataset['Name'], 
     'VCF [File]'=>File.join(@result_dir, "#{@dataset['Name']}.vcf"), 
     'Species'=>@dataset['Species'],
     'build'=>@params['build'],
     'paired'=>@params['paired'],
     'Read Count'=>@dataset['Read Count']
    }.merge factor_dataset
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
    command << "input = list()\n"
    input = @dataset
    input.keys.each do |key|
      command << "input[['#{key}']] = '#{input[key]}'\n" 
    end
    command << "output = list()\n"
    output = next_dataset
    output.keys.each do |key|
      command << "output[['#{key}']] = '#{output[key]}'\n" 
    end
    command << "pileup2vcf(input=input, output=output, config=config)\n"
    command << "EOT"
    command
  end
end

if __FILE__ == $0
  run PileupApp
  #usecase = Bowtie2App.new

  #usecase.project = "p1001"
  #usecase.user = 'masamasa'

  # set user parameter
  # for GUI sushi
  #usecase.params['process_mode'].value = 'SAMPLE'
  #usecase.params['build'] = 'mm10'
  #usecase.params['paired'] = true
  #usecase.params['strandMode'] = 'both'
  #usecase.params['cores'] = 8
  #usecase.params['node'] = 'fgcz-c-048'

  # also possible to load a parameterset csv file
  # mainly for CUI sushi
  #usecase.parameterset_tsv_file = 'tophat_parameterset.tsv'
  #usecase.parameterset_tsv_file = 'test.tsv'

  # set input dataset
  # mainly for CUI sushi
  #usecase.dataset_tsv_file = 'tophat_dataset.tsv'

  # also possible to load a input dataset from Sushi DB
  #usecase.dataset_sushi_id = 3

  # run (submit to workflow_manager)
  #usecase.run
  #usecase.test_run

end

