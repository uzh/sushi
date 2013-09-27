#!/usr/bin/env ruby
# encoding: utf-8

require 'sushiApp'

class ProtvizApp < SushiApp
  def initialize
    super
    @name = 'Protviz'
    @params['process_mode'] = 'SAMPLE'
    @analysis_category = 'Demo'
    @required_columns = ['Name','Mascot DAT']
    @required_params = ['name']
    @params['cores'] = '1'
    @params['ram'] = '2'
    @params['scratch'] = '10'
    @params['name'] = 'protViz_MarkerFinder'
    @params['markerIonNames'] = ['manual', 'ADP_Ribose', 'Glycan', 'HexNAc' ]
    @params['markerIonMasses'] = ''
    @params['itol_ppm'] = 20
  end
  def next_dataset
    {'Name'=>@params['name'],
     'MGF [File]'=>File.join(@result_dir, "#{@dataset['Name']}.mgf"),
     'PDF [File,Link]'=>File.join(@result_dir, "#{@dataset['Name']}.pdf"),
     'CSV [File]'=>File.join(@result_dir, "#{@dataset['Name']}.csv")
    }
  end
  def commands
    command =<<-EOS
/srv/GT/analysis/hubert/mascotProtviz/mascotDat2RData.pl -d=#{File.join(@gstore_dir, @dataset['Mascot DAT'])} -m=/srv/GT/analysis/hubert/mascotProtviz/mod_file

/usr/local/ngseq/bin/R --vanilla --slave << EOT
dataFileBase = sub(".dat\$", "", basename("#{@dataset['Mascot DAT']}"))
markerIonsName = "#{@params['markerIonNames']}"
itol_ppm = #{@params['itol_ppm']}
mgfFile = "#{@dataset['Name']}.mgf"
pdfFile = "#{@dataset['Name']}.pdf"
csvFile = "#{@dataset['Name']}.csv"
ionMasses = c(#{@params['markerIonMasses']})

library(protViz)
markerIons = switch(markerIonsName,
                    Glycan=sort(c(366.140, 325.1135, 274.0927, 204.0866, 163.0601)),
                    ADP_Ribose=sort(c(428.0367, 348.0704, 250.0935, 136.0618)),
                    HexNAc=sort(c(204.0866, 186.0761, 168.0655, 144.0655, 138.0550, 126.0550)),
                    manual=sort(ionMasses))
PTM_MarkerFinder_util(dataFileBase, 
    mZmarkerIons=markerIons, 
    minMarkerIntensityRatio=2,
    minNumberIons=2, 
    itol_ppm=itol_ppm, 
    write_csv=TRUE)
file.rename(paste(dataFileBase, "mgf", sep="."), mgfFile)
file.rename(paste(dataFileBase, "pdf", sep="."), pdfFile)
file.rename(paste(dataFileBase, "csv", sep="."), csvFile)
EOT
EOS
  end
end

if __FILE__ == $0
  usecase = ProtvizApp.new

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

