#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'

class KmeansDemoSampleApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Kmeans per Sample Demo'
    @params['process_mode'] = 'SAMPLE'
    @analysis_category = 'Clustering'
    @required_columns = ['Name','Count']
    @required_params = ['name', 'nClusters']
    @params['cores'] = '1'
    @params['ram'] = '2'
    @params['scratch'] = '10'
    @params['name'] = 'k-means_clustering'
    @params['nClusters'] = '5'
    @params['makePlot'] = false
  end
  def next_dataset
    {'Name'=>@params['name'],
     'Cluster [File]'=>File.join(@result_dir, "#{@dataset['Name']}.txt")
    }
  end
  def commands
    command =<<-EOS
/usr/local/ngseq/bin/R --vanilla --slave << EOT
nClusters = #{@params['nClusters']}
countFile = "#{File.join(@gstore_dir, @dataset['Count'])}"
clusterFile = basename("#{next_dataset['Cluster [File]']}")
makePlot = "#{@params['makePlot']}" == "true"
x = read.table(countFile, sep="\\t", stringsAsFactors=FALSE, header=TRUE)
logCounts = log2(x[ , "transcriptCount"] + 1)
kmeansResult = kmeans(logCounts, nClusters)
clusterTable = data.frame(Name=rownames(x), Cluster=kmeansResult[["cluster"]])
write.table(clusterTable, file=clusterFile, row.names=FALSE)
#if (makePlot){
#    png(file=basename(output$"CenterPlot [File]"), height=800, width=800)
#    my.profilePlot(cl$centers)
#    dev.off()
#  }
EOT
EOS
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

