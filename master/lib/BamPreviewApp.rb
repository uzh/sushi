#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class BamPreviewApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'BAM Preview'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'QC'
@description =<<-EOS
    Run a mapper and compute stats on the bam files<br/>
    If no <i>mapOptions</i> are provided the default map options are:<ul>
<li> STAR: --outFilterType BySJout --outFilterMatchNmin 30 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.05 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignIntronMax 1000000 --alignMatesGapMax 1000000  --outFilterMultimapNmax 50 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --outSAMstrandField intronMotif</li>
<li> bowtie: ''</li>
<li> bowtie2: --no-unal </li>
<li> tophat: --mate-inner-dist 100 --mate-std-dev 150 </li>
<li> bwa-mem: ''</li>
</ul>
EOS
    @required_columns = ['Name','Read1','Species']
    @required_params = ['name', 'paired']
    @params['cores'] = '8'
    @params['ram'] = '50'
    @params['scratch'] = '100'
    @params['paired'] = false
    @params['name'] = 'BAM_Preview'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['subsampleReads'] = 30
    @params['subsampleReads', 'description'] = 'take only every nth read'
    @params['mapMethod'] = ['STAR', 'bowtie', 'bowtie2', 'tophat', 'bwa-mem']
    @params['mapOptions'] = ''
    @params['trimAdapter'] = true
    @params['trimLeft'] = 0
    @params['trimRight'] = 0
    @params['minTailQuality'] = 0
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["QC/Trimmomatic", "QC/Flexbar", "Tools/samtools", "Aligner/Bowtie2", "Aligner/STAR", "Dev/Python2", "Dev/R", "Tools/sambamba"]
  end
  def next_dataset
    report_dir = File.join(@result_dir, @params['name'])
    {'Name'=>@params['name'],
     'Report [File]'=>report_dir,
     'Html [Link]'=>File.join(report_dir, '00index.html'),
     'Species'=>(dataset = @dataset.first and dataset['Species']),
     'refBuild'=>@params['refBuild'],
     'refFeatureFile'=>@params['refFeatureFile']
    }
  end
  def set_default_parameters
    @params['paired'] = dataset_has_column?('Read2')
    if dataset_has_column?('strandMode')
      @params['strandMode'] = @dataset[0]['strandMode']
    end
  end

  def commands
    run_RApp("EzAppBamPreview")
  end
end

if __FILE__ == $0
  usecase = BamStatsApp.new

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
