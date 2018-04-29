#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-100301'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class VirDetectApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'VirDetect'
    @analysis_category = 'Metagenomics'
    @description =<<-EOS
Virome analysis for the diagnostic use in veterinary virology<br/>
<a href='http://www.vetvir.uzh.ch/en/Research/Virology/Molecular-and-Clinical-Veterinary-Virology/Virome-Analysis.html'>Virome Analysis Group</a>, The Institute of Virology, UZH<br>
<a href='https://www.ncbi.nlm.nih.gov/pubmed/25009045'>Reference paper</a>, base on which this pileline is implemented.
EOS

    @required_columns = ['Name','Read1','Species']
    @required_params = ['virBuild','hostBuild','paired']
    @params['virBuild'] = ref_selector
    @params['virBuild', 'description'] = 'the viral reference database to use.'
    @params['hostBuild'] = ref_selector
    @params['hostBuild', 'description'] = 'the none-human host genome to use as reference for removal of contamination'
    @params['paired'] = false
    @params['paired', 'description'] = 'whether the reads are paired end; if false then only Read1 is considered even if Read2 is available.'
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '60'
    @params['scratch'] = '150'
    @params['cmdOptionsHost'] = '--very-sensitive'
    @params['cmdOptionsHost', 'description'] = 'specify the commandline options for bowtie2 during mapping to contaminated genome(s); do not specify any option that is already covered by the dedicated input fields'
    @params['cmdOptions'] = '-a --very-sensitive --no-mixed --no-discordant -X 1000'
    @params['cmdOptions', 'description'] = 'specify the commandline options for bowtie2 during mapping to the viral reference database; do not specify any option that is already covered by the dedicated input fields'
    @params['trimAdapter'] = true
    @params['trimAdapter', 'description'] = 'if adapters should be trimmed'
    @params['trimLeft'] = 5
    @params['trimLeft', 'description'] = 'fixed trimming at the "left" i.e. 5-prime end of the read'
    @params['trimRight'] = 0
    @params['trimRight', 'description'] = 'fixed trimming at the "right" i.e. 3-prime end of the read'
    @params['minTailQuality'] = 10
    @params['minTailQuality', 'description'] = 'if above zero, then reads are trimmed as soon as 4 consecutive bases have lower mean quality'
    @params['minAvgQuality'] = 20
    @params['minReadLength'] = 50
    @params['minReadCount'] = 10
    @params['minReadCount', 'description'] = 'minimum number of reads mapped to a given genome so that it is reported as detected.'
    @params['specialOptions'] = ''
    @params['specialOptions', 'description'] = 'special unsupported options that the R wrapper may support, format: <key>=<value>'
    @params['mail'] = ""
    @modules = ["Tools/samtools", "Aligner/Bowtie2", "QC/Flexbar", "QC/Trimmomatic", "Tools/BEDTools", "Dev/R", "Tools/sambamba"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
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
     'Species'=>@dataset['Species'],
     'virBuild'=>@params['virBuild'],
     'hostBuild'=>@params['hostBuild'],
     'paired'=>@params['paired'],
     'Read Count'=>@dataset['Read Count'],
     'OutDir [File]'=>File.join(@result_dir, "#{@dataset['Name']}"),
     'OutReport [Link]'=>File.join(@result_dir, "#{@dataset['Name']}", "#{@dataset['Name']}.html")
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppVirDetect")
  end
end

if __FILE__ == $0
  run Bowtie2App
  #usecase = VirDetectApp.new

  #usecase.project = "p1001"
  #usecase.user = 'masamasa'

  # set user parameter
  # for GUI sushi
  #usecase.params['process_mode'].value = 'SAMPLE'
  #usecase.params['refBuild'] = 'mm10'
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
