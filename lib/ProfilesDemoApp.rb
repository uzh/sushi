#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'

class ProfilesDemoApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Profile Plot'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Visualization'
    @required_columns = ['Name','Count']
    @required_params = ['name', 'doLogTransform']
    @params['cores'] = '1'
    @params['ram'] = '2'
    @params['scratch'] = '10'
    @params['name'] = 'profile_plot'
    @params['doLogTransform'] = true
  end
  def next_dataset
    {'Name'=>@params['name'],
     'Profiles [File,Link]'=>File.join(@result_dir, "profiles.png")
    }
  end
  def commands
    command =<<-EOS
 R --vanilla --slave << EOT
 doLogTransform = "#{@params['doLogTransform']}" == "true"
 pngFile =  basename("#{next_dataset['Profiles [File,Link]']}")
 inputDatasetFile = '#{@input_dataset_tsv_path}'
 dataset = read.table(inputDatasetFile, sep="\\t", stringsAsFactors=FALSE, check.names=FALSE, header=TRUE)
 countFiles = file.path('#{@gstore_dir}', dataset[["Count [File]"]])
 countList = lapply(countFiles, function(cf){ read.table(cf, sep="\\t", stringsAsFactors=FALSE, check.names=FALSE, header=TRUE)[["transcriptCount"]]})
 if (doLogTransform){
  countList = lapply(countList, function(x){log(x+1)})
  bw = 0.1
 } else {
  bw = 100
 }
 densList = lapply(countList, density, bw=bw)
 yMax = max(sapply(densList, function(x){max(x[['y']])}))
 xMax = max(sapply(densList, function(x){max(x[['x']])}))
 png(file=pngFile, height=800, width=800) 
 plot(1, type="n", xlim=c(0.5, xMax), ylim=c(0, yMax), xlab="counts", ylab="freq")
 sapply(densList, function(x){lines(x[['x']], x[['y']])})
 dev.off()
EOT
EOS
  end
end

if __FILE__ == $0
  usecase = ProfilesDemoApp.new

  usecase.project = "p1001"
  usecase.user = 'masamasa'

  usecase.params['name'] = "profiles_testrun"
  usecase.params['doLogTransform'] = true

  # also possible to load a parameterset csv file
  # mainly for CUI sushi
  #usecase.parameterset_tsv_file = 'tophat_parameterset.tsv'

  # set input dataset
  # mainly for CUI sushi
  usecase.dataset_tsv_file = '/Users/hubert/sushi/Demo-readcounts.tsv'

  # also possible to load a input dataset from Sushi DB
  #usecase.dataset_sushi_id = 3

  # run (submit to workflow_manager)
  #usecase.run
  usecase.test_run

end

