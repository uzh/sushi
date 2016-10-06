#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class GatkDnaHaplotyperApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'GatkDnaHaplotyper'
    @analysis_category = 'Variants'
    @description =<<-EOS
Haplotype calling for DNA-seq<br/>
<a href='https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php'>HaplotypeCaller</a>
    EOS
    @required_columns = ['Name','BAM','BAI', 'refBuild']
    @required_params = ['name']
    @params['cores'] = '8'
    @params['ram'] = '50'
    @params['scratch'] = '100'
    @params['name'] = 'GATK_DnaVariants'
    @params['refBuild'] = ref_selector
    @params['targetFile'] = '/srv/GT/analysis/lopitz/GATK-tutorial_data/intervals/test.bed'
    @params['getRealignedBam'] = true
    @params['markDuplicates'] = false
    @params['addReadGroup'] = false
    @params['specialOptions'] = ''
    @params['mail'] = ""
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'GVCF [File]'=>File.join(@result_dir, "#{@dataset['Name']}-HC_calls.g.vcf.gz"),
     'BAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}-realigned.bam"),
     'BAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}-realigned.bai"),
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild']
    }.merge(extract_column("Factor")).merge(extract_column("B-Fabric"))
  end
  def commands
    run_RApp("EzAppGatkDnaHaplotyper")
  end
end

if __FILE__ == $0

end
