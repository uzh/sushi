#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class GatkRnaHaplotyperApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'GatkRnaHaplotyper'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Variants'
    @description =<<-EOS
Haplotype calling for RNA-seq<br/>
<a href='https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php'>HaplotypeCaller</a>
    EOS
    @required_columns = ['Name','BAM','BAI', 'build', 'Species']
    @required_params = ['name', 'paired']
    @params['cores'] = '24'
    @params['ram'] = '100'
    @params['scratch'] = '500'
    @params['paired'] = false
    @params['name'] = 'GATK_RnaVariants'
    @params['build'] = ref_selector
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Tools/samtools", "Dev/jdk", "Variants/GATK", "Tools/Picard", "Dev/R", "Tools/sambamba"]
  end
  def next_dataset
    report_dir = File.join(@result_dir, @params['name'])
    {'Name'=>@params['name'],
     'VCF [File]'=>File.join(@result_dir, "#{@params['name']}-haplo.vcf.gz"),
     'TBI [File]'=>File.join(@result_dir, "#{@params['name']}-haplo.vcf.gz.tbi"),
     'Report [File]'=>report_dir,
     'Html [Link]'=>File.join(report_dir, '00index.html'),
     'Species'=>(dataset = @dataset.first and dataset['Species']),
     'build'=>@params['build']
    }
  end
  def set_default_parameters
    @params['build'] = @dataset[0]['build']
    if dataset_has_column?('paired')
      @params['paired'] = @dataset[0]['paired']
    end
  end

  def commands
    run_RApp('EzAppGatkRnaHaplotyper')
  end
end

if __FILE__ == $0
  usecase = BamStatsApp.new

  usecase.project = "p1001"
  usecase.user = "masa"

  # set user parameter
  # for GUI sushi
  #usecase.params['process_mode'].value = 'SAMPLE'
  #usecase.params['build'] = 'TAIR10'
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
