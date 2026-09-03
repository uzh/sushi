#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class MpileupApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'samtools mpileup'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Variants'
    @description =<<-EOS
Variant analysis with samtools/bcftools.<br/>
The analysis runs the 3 commands of bcftools: mpileup, call, and filter<br/>
<a href='http://www.htslib.org/doc/bcftools.html'>bcftools manual/</a><br>
EOS
    @required_columns = ['Name','BAM','BAI', 'refBuild', 'Species']
    @required_params = ['name', 'paired']
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '30'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '100'
    @params['scratch', "context"] = "slurm"
    @params['paired'] = false
    @params['paired', "context"] = "Mpileup"
    @params['name'] = 'Mpileup_Variants'
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "referfence genome assembly"
    @params['region'] = ""
    @params['region', 'description'] = 'The region of the genome. You can give either a chromosome name or a region on a chromosome like chr1:1000-2000'
    @params['region', "context"] = "Mpileup"
    @params['mpileupOptions'] = '--skip-indels --annotate AD,INFO/AD,ADF,ADR,SP'
    @params['mpileupOptions', 'description'] = 'The options to the bcftools mpileup command'
    @params['mpileupOptions', "context"] = "Mpileup"
    @params['callOptions'] = '--multiallelic-caller --keep-alts --variants-only'
    @params['callOptions', 'description'] = 'The options to <a href=http://www.htslib.org/doc/bcftools.html#call>bcftools call</a>'
    @params['callOptions', "context"] = "Mpileup"
    @params['filterOptions'] = '--include "MIN(DP)>5"'
    @params['filterOptions', 'description'] = 'The options to <a href=http://www.htslib.org/doc/bcftools.html#filter>bcftools filter</a>'
    @params['filterOptions', "context"] = "Mpileup"
    @params['specialOptions'] = ''
    @params['specialOptions', 'description'] = 'special unsupported options that the R wrapper may support, format: <key>=<value>'
    @params['mail'] = ""
    @modules = ["Tools/samtools", "Tools/bcftools", "Dev/jdk", "Tools/Picard", "Dev/R"]
    @inherit_columns = ["Order Id"]
  end
  def next_dataset
    report_dir = File.join(@result_dir, @params['name'])
    {'Name'=>@params['name'],
     'VCF [File]'=>File.join(@result_dir, "#{@params['name']}.vcf.gz"),
     'TBI [File]'=>File.join(@result_dir, "#{@params['name']}.vcf.gz.tbi"),
     #'IGV Starter [Link]'=>File.join(@result_dir, "#{@params['name']}-igv.jnlp"),
     'Report [File]'=>report_dir,
     'Html [Link]'=>File.join(report_dir, '00index.html'),
     'Species'=>(dataset = @dataset.first and dataset['Species']),
     'refBuild'=>@params['refBuild']
     #'IGV Starter [File]'=>File.join(@result_dir, "#{@params['name']}-igv.jnlp"),
     #'IGV Session [File]'=>File.join(@result_dir, "#{@params['name']}-igv.xml")
    }.merge(extract_columns(colnames: @inherit_columns))
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('paired')
      @params['paired'] = @dataset[0]['paired']
    end
  end

  def commands
    run_RApp("EzAppMpileup")
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
