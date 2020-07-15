#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class TeqcApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'Teqc'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'QC'
    @description =<<-EOS
Quality control for target capture experiments<br/>
<a href='http://www.bioconductor.org/packages/release/bioc/html/TEQC.html'>manual</a>
    EOS
    @required_columns = ['Name','BAM','BAI', 'refBuild']
    @required_params = ['name', 'paired']
    @params['cores'] = '4'
    @params['ram'] = '100'
    @params['scratch'] = '100'
    @params['paired'] = false
    @params['name'] = 'TEQC_Result'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['designFile'] = {'select'=>''}
    Dir["/srv/GT/databases/targetEnrichment_designs/*"].sort.select{|design| File.directory?(design)}.each do |dir|
      @params['designFile'][File.basename(dir)] = File.basename(dir)
    end
    @params['removeDuplicates'] = true
    @params['covUniformityPlot'] = true
    @params['covTargetLengthPlot'] = true
    @params['duplicatesPlot'] = true
    @params['cmdOptions'] = ""
    @params['mail'] = ""
    @modules = ["Tools/samtools", "Tools/sambamba","Dev/R"]
  end
 def set_default_parameters
   @params['refBuild'] = @dataset[0]['refBuild']
   if dataset_has_column?('refFeatureFile')
     @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
   end
   if dataset_has_column?('paired')
      @params['paired'] = @dataset[0]['paired']
    end
 end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@params['name'],
     'Report [File]'=>report_file,
     'Html [Link]'=>report_link,
    }
  end
  def commands
    run_RApp("EzAppTeqc")
  end
end


if __FILE__ == $0
  usecase = TeqcApp.new

  usecase.project = "p1001"
  usecase.user = "masa"

  # set user parameter
  # for GUI sushi
  #usecase.params['process_mode'].value = 'SAMPLE'
  #usecase.params['refBuild'] = 'TAIR10'
  #usecase.params['paired'] = true
  #usecase.params['cores'] = 2
  #usecase.params['node'] = 'fgcz-c-048'

  # also possible to load a parameterset csv file
  # mainly for CUI sushi
  usecase.parameterset_tsv_file = 'tophat_parameterset.tsv'
  #usecase.params['name'] = 'name'

  # set input dataset
  # mainly for CUI sushi
  usecase.dataset_tsv_file = 'tophat_dataset.tsv'

  # also possible to load a input dataset from Sushi DB
  #usecase.dataset_sushi_id = 1

  # run (submit to workflow_manager)
  usecase.run
  #usecase.test_run

end
