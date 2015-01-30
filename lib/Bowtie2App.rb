#!/usr/bin/env ruby
# encoding: utf-8
Version = '20150126-164943'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class Bowtie2App < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Bowtie2'
    @analysis_category = 'Map'
    @description =<<-EOS
Fast and sensitive read alignment. Supports local and end-to-end mode<br/>
<a href='http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml'>http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml/</a>
EOS
    
    @required_columns = ['Name','Read1','Species']
    @required_params = ['build','paired']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '16'
    @params['scratch'] = '100'
    @params['build'] = ref_selector
    @params['build', 'description'] = 'the genome build and annotation to use as reference. If human variant calling is the main goal, please use hg_19_karyotypic.'
    @params['paired'] = false
    @params['paired', 'description'] = 'whether the reads are paired end; if false then only Read1 is considered even if Read2 is available.'
    @params['cmdOptions'] = '--no-unal'
    @params['cmdOptions', 'description'] = 'specify the commandline options for bowtie2; do not specify any option that is already covered by the dedicated input fields'
    @params['trimAdapter'] = false
    @params['trimAdapter', 'description'] = 'if adapters should be trimmed'
    @params['trimLeft'] = 0
    @params['trimLeft', 'description'] = 'fixed trimming at the "left" i.e. 5-prime end of the read'
    @params['trimRight'] = 0
    @params['trimRight', 'description'] = 'fixed trimming at the "right" i.e. 3-prime end of the read'
    @params['minTailQuality'] = 0
    @params['minTailQuality', 'description'] = 'if above zero, then reads are trimmed as soon as 4 consecutive bases have lower mean quality'
    #@params['specialOptions'] = ''
    #@params['specialOptions', 'description'] = 'special unsupported options that the R wrapper may support, format: <key>=<value>'
    @params['mail'] = ""
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
    {'Name'=>@dataset['Name'], 
     'BAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam"), 
     'BAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam.bai"),
     'IGV Starter [Link]'=>File.join(@result_dir, "#{@dataset['Name']}-igv.jnlp"),
     'Species'=>@dataset['Species'],
     'build'=>@params['build'],
     'paired'=>@params['paired'],
     'Read Count'=>@dataset['Read Count'],
     'IGV Starter [File]'=>File.join(@result_dir, "#{@dataset['Name']}-igv.jnlp"),
     'IGV Session [File]'=>File.join(@result_dir, "#{@dataset['Name']}-igv.xml")
    }.merge(extract_column("Factor")).merge(extract_column("B-Fabric"))
  end
  def commands
    command = "/usr/local/ngseq/bin/R --vanilla --slave << EOT\n"
    command << "GLOBAL_VARIABLES <<- '#{GLOBAL_VARIABLES}'\n"
    command << "R_SCRIPT_DIR <<- '#{R_SCRIPT_DIR}'\n"
    command<<  "source(file.path(R_SCRIPT_DIR, 'init.R'))\n"
    command << "config = list()\n"
    config = @params
    config.keys.each do |key|
      command << "config[['#{key}']] = '#{config[key]}'\n" 
    end
    command << "config[['dataRoot']] = '#{@gstore_dir}'\n"
    command << "config[['resultDir']] = '#{@result_dir}'\n"
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
    command << "runApp('mapBowtie2App', input=input, output=output, config=config)\n"
    command << "EOT"
    command
  end
end

if __FILE__ == $0
  run Bowtie2App
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

