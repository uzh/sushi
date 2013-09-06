#!/usr/bin/env ruby
# encoding: utf-8

require 'sushiApp'

class KmeansApp < SushiApp
  def initialize
    super
    @name = 'Kmeans'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Clustering'
    @required_columns = ['Name','Count']
    @required_params = ['name', 'nClusters']
    @params['cores'] = '1'
    @params['ram'] = '2'
    @params['scratch'] = '10'
    @params['name'] = 'k-means_clustering'
    @params['nClusters'] = '5'
  end
  def next_dataset
    {'Name'=>@params['name'],
     'Cluster [File]'=>File.join(@result_dir, "#{@params['name']}.txt")
    }
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
    command << "kmeansDemoApp(input=inputDatasetFile, output=output, config=config)\n"
    command << "EOT"
    command
  end
end

if __FILE__ == $0
  usecase = KmeansApp.new

  usecase.project = "p1001"
  usecase.user = 'masamasa'

  usecase.params['name'] = "kmeans_testrun"
  usecase.params['nClusters'] = 5

  # also possible to load a parameterset csv file
  # mainly for CUI sushi
  #usecase.parameterset_tsv_file = 'tophat_parameterset.tsv'

  # set input dataset
  # mainly for CUI sushi
  usecase.dataset_tsv_file = '/srv/GT/analysis/hubert/sushi/p1001-countDataset.txt'

  # also possible to load a input dataset from Sushi DB
  #usecase.dataset_sushi_id = 3

  # run (submit to workflow_manager)
  #usecase.run
  usecase.test_run

end

