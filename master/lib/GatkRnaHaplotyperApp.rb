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
    @required_columns = ['Name','BAM','BAI', 'refBuild', 'Species']
    @required_params = ['name', 'paired']
    @params['cores'] = '4'
    @params['ram'] = '30'
    @params['scratch'] = '500'
    @params['paired'] = false
    @params['name'] = 'GATK_RnaVariants'
    @params['refBuild'] = ref_selector
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Tools/samtools", "Dev/jdk", "Variants/GATK", "Tools/Picard", "Dev/R"]
  end
  def next_dataset
    report_dir = File.join(@result_dir, @params['name'])
    {'Name'=>@params['name'],
     'VCF [File]'=>File.join(@result_dir, "#{@params['name']}-haplo.vcf.gz"),
     'TBI [File]'=>File.join(@result_dir, "#{@params['name']}-haplo.vcf.gz.tbi"),
     'Report [File]'=>report_dir,
     'Html [Link]'=>File.join(report_dir, '00index.html'),
     'Species'=>(dataset = @dataset.first and dataset['Species']),
     'refBuild'=>@params['refBuild']
    }
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('paired')
      @params['refBuild'] = @dataset[0]['refBuild']
    end
  end

  def commands
    run_RApp('EzAppGatkRnaHaplotyper')
  end
end

