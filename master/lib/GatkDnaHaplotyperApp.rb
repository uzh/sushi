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
    @params['targetFile'] = {'select'=>''}
    if defined?(TARGET_ENRICHMENT_DESIGN_DIR)
      Dir["#{TARGET_ENRICHMENT_DESIGN_DIR}/*.bed"].sort.select{|bed| File.file?(bed)}.each do |file|
        @params['targetFile'][File.basename(file)] = File.basename(file)
      end
    end
    @params['getRealignedBam'] = false
    @params['markDuplicates'] = false
    @params['addReadGroup'] = true
    @params['knownSitesAvailable'] = false
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Tools/samtools", "Dev/jdk", "Variants/GATK", "Tools/Picard", "Tools/htslib", "Dev/R", "Tools/sambamba"]
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
      'targetFile'=>@params['targetFile'],
      'refBuild'=>@params['refBuild']
    }.merge(extract_columns(@inherit_tags))

    if @params['getRealignedBam']
      dataset['BAM [File]'] = File.join(@result_dir, "#{@dataset['Name']}-realigned.bam")
      dataset['BAI [File]'] = File.join(@result_dir, "#{@dataset['Name']}-realigned.bai")
    end
    dataset
  end
  def commands
    run_RApp("EzAppGatkDnaHaplotyper")
  end
end

if __FILE__ == $0

end
