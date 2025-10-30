#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class Mutect2App <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'Mutect2'
    @analysis_category = 'Variants'
    @description =<<-EOS
Somatic variant calling for DNA-seq<br/>
<a href='https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2'>Mutect2</a>
    EOS
    @required_columns = ['Name','BAM','BAI', 'CtrlBam', 'refBuild']
    @required_params = ['name']
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '50'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '100'
    @params['scratch', "context"] = "slurm"
    @params['name'] = 'GATK_somaticDnaVariants'
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "referfence genome assembly"
    @params['snpEffDB'] = ''
    @params['snpEffDB', "context"] = "Mutect2"
    @params['TumorOnlyMode'] = false
    @params['TumorOnlyMode', 'description'] = 'default mode is tumor with matched normal,  TumorOnlyMode=true ignores the normal sample'
    @params['TumorOnlyMode', "context"] = "Mutect2"
    @params['specialOptions'] = ''
    @params['cmdOptions'] = ''
    @params['cmdOptions', "context"] = "Mutect2"
    @params['mail'] = ""
    @modules = ["Dev/jdk", "Variants/GATK", "Variants/SnpEff/4.3","Dev/R", "Tools/Picard", "Tools/samtools"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
  end
  def next_dataset
    dataset = {
    'Name'=>@dataset['Name'],
      'VCF [File]'=>File.join(@result_dir, "#{@dataset['Name']}.somatic.ann.vcf.gz"),
      'VCFINDEX [File]'=>File.join(@result_dir, "#{@dataset['Name']}.somatic.ann.vcf.gz.tbi"),
      'Other [File]'=>File.join(@result_dir, "#{@dataset['Name']}_misc.zip"),
      'Species'=>@dataset['Species'],
      'refBuild'=>@params['refBuild']
    }.merge(extract_columns(@inherit_tags))
    dataset
  end
  def commands
    run_RApp("EzAppMutect2")
  end
end

if __FILE__ == $0

end
