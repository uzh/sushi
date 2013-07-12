#!/usr/bin/env ruby
# encoding: utf-8
Version = '20130712-065915'

require 'sushiApp'

class Bowtie2WrappedApp < SushiApp
  def initialize
    super
    @name = 'Bowtie2 Wrapped'
    @analysis_category = 'Map'
    @required_columns = ['Sample','Read1','Species']
    @required_params = ['build','paired', 'strandMode']
    # optional params
    @params['build'] = {'select'=>''}
    Dir["/srv/GT/reference/*/*/*"].sort.select{|build| File.directory?(build)}.each do |dir|
      @params['build'][dir.gsub(/\/srv\/GT\/reference\//,'')] = File.basename(dir)
    end
    @params['paired'] = false
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['featureFile'] = 'genes.gtf'
    @params['cmdOptions'] = ''
    @params['trimAdapter'] = false
    @params['trimLeft'] = 0
    @params['trimRight'] = 0
    @params['minTailQuality'] = 0
    @params['specialOptions'] = ''
    @output_files = ['BAM','BAI']
  end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  end
  def next_dataset
    {'Sample'=>@dataset['Sample'], 
     'BAM'=>File.join(@result_dir, "#{@dataset['Sample']}.bam"), 
     'BAI'=>File.join(@result_dir, "#{@dataset['Sample']}.bam.bai"),
     'Build'=>@params['build']
    }
  end
  def commands
    command = "/usr/local/ngseq/bin/R --vanilla --slave << EOT\n"
    command << "source('/usr/local/ngseq/sushi_scripts/init.R')\n"
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
    command << "mapBowtie2(input=input, output=output, config=config)\n"
    command << "EOT"
    command
  end
end

if __FILE__ == $0
  usecase = Bowtie2WrappedApp.new

  usecase.project = "p1001"
  usecase.user = 'masamasa'

  # set user parameter
  # for GUI sushi
  #usecase.params['process_mode'].value = 'SAMPLE'
  usecase.params['build'] = 'mm10'
  usecase.params['paired'] = true
  usecase.params['strandMode'] = 'both'
  usecase.params['cores'] = 8
  usecase.params['node'] = 'fgcz-c-048'

  # also possible to load a parameterset csv file
  # mainly for CUI sushi
  #usecase.parameterset_tsv_file = 'tophat_parameterset.tsv'
  #usecase.parameterset_tsv_file = 'test.tsv'

  # set input dataset
  # mainly for CUI sushi
  #usecase.dataset_tsv_file = 'tophat_dataset.tsv'

  # also possible to load a input dataset from Sushi DB
  usecase.dataset_sushi_id = 3

  # run (submit to workflow_manager)
  usecase.run
  #usecase.test_run

end

