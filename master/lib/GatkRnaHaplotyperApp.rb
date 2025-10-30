#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class GatkRnaHaplotyperApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'GatkRnaHaplotyper'
    @analysis_category = 'Variants'
    @description =<<-EOS
Haplotype calling for RNA-seq<br/>
<a href='https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-'>RNA-seq Best practice</a>
    EOS
    @required_columns = ['Name','BAM','BAI', 'refBuild']
    @required_params = ['name']
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '50'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '100'
    @params['scratch', "context"] = "slurm"
    @params['name'] = 'GATK_RnaVariants'
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "referfence genome assembly"
    @params['refFeatureFile'] = 'genes.gtf'
    @params['refFeatureFile', "context"] = "GatkRnaHaplotyper"
    #@params['getRealignedBam'] = false
    @params['markDuplicates'] = false
    @params['markDuplicates', "context"] = "Picard"
    @params['addReadGroup'] = true
    @params['addReadGroup', "context"] = "GatkRnaHaplotyper"
    @params['dbsnpFile'] = ''
    @params['dbsnpFile', "context"] = "GatkRnaHaplotyper"
    @params['specialOptions'] = ''
    @params['specialOptions', "context"] = "GatkRnaHaplotyper"
    @params['mail'] = ""
    @modules = ["Tools/samtools", "Variants/GATK", "Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
  end
  def next_dataset
    dataset = {
    'Name'=>@dataset['Name'],
      'GVCF [File]'=>File.join(@result_dir, "#{@dataset['Name']}-HC_calls.g.vcf.gz"),
      'GVCFINDEX [File]'=>File.join(@result_dir, "#{@dataset['Name']}-HC_calls.g.vcf.gz.tbi"),
      'Species'=>@dataset['Species'],
      'refBuild'=>@params['refBuild']
    }.merge(extract_columns(@inherit_tags))

#    if @params['getRealignedBam']
#      dataset['BAM [File]'] = File.join(@result_dir, "#{@dataset['Name']}-realigned.bam")
#      dataset['BAI [File]'] = File.join(@result_dir, "#{@dataset['Name']}-realigned.bai")
#    end
    dataset
  end
  def commands
    run_RApp("EzAppGatkRnaHaplotyper")
  end
end

if __FILE__ == $0

end
