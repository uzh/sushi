#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class GATKv4DNAHaplotypeCallerGVCFApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'GATKv4DNAHaplotypeCallerGVCF'
    @analysis_category = 'Variants'
    @description =<<-EOS
Haplotype calling for DNA-seq with > version 4.0 in GVCF mode<br/>
<a href='https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.8.0/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php'>HaplotypeCaller</a>
    EOS
    @required_columns = ['Name', 'BAM', 'refBuild', "Dummy"]
    @required_params = ['name']
    @params['cores'] = '4'
    @params['ram'] = '50'
    @params['scratch'] = '100'
    @params['name'] = 'GATKv4_gVCF'
    @params['refBuild'] = ref_selector
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Tools/samtools/1.9", "Variants/GATK/4.1.2.0", "Tools/Picard/2.18.0"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic", "BAM"]
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
  end
  def next_dataset
    dataset = {
      'Name'=>@dataset['Name'],
      'GVCF [File]'=>File.join(@result_dir, "#{@dataset['Name']}.g.vcf.gz"),
      'GVCFINDEX [File]'=>File.join(@result_dir, "#{@dataset['Name']}.g.vcf.idx"),
      'Species'=>@dataset['Species'],
      'refBuild'=>@params['refBuild'],
      'Dummy [File]'=>File.join(@result_dir, "#{@dataset['Name']}_dummy.txt")
    }.merge(extract_columns(@inherit_tags))

    dataset
  end
  def commands
    ref = File.join(GENOME_REF_DIR, @dataset['refBuild'].split('/')[0,3].join("/")+"/Sequence/WholeGenomeFasta/genome.fa")
    bam = File.join(@gstore_dir, @dataset['BAM'])
    sort_bam = File.basename(@dataset['BAM']).gsub(/.bam/, '.sort.bam') 
    sort_rg_bam = File.basename(@dataset['BAM']).gsub(/.bam/, '.sort.rg.bam')
    sort_rg_dup_bam = File.basename(@dataset['BAM']).gsub(/.bam/, '.sort.rg.dup.bam')
    sort_rg_dup_met = File.basename(@dataset['BAM']).gsub(/.bam/, '.sort.rg.dup.met')
    g_vcf = @dataset['Name'] + '.g.vcf'
    command = "samtools sort #{bam} -o #{sort_bam}\n"
    command << "gatk AddOrReplaceReadGroups -I #{sort_bam} -ID #{@dataset['Name']} -PU #{@dataset['Name']} -LB #{@dataset['Name']} -SM #{@dataset['Name']} -PL Illumina -O #{sort_rg_bam}\n"
    command << "samtools index #{sort_rg_bam}\n"
    command << "gatk MarkDuplicates -I #{sort_rg_bam} -O #{sort_rg_dup_bam} -M #{sort_rg_dup_met}\n"
    command << "samtools index #{sort_rg_dup_bam}\n"
    command << "gatk HaplotypeCaller -R #{ref} -I #{sort_rg_dup_bam} -O #{g_vcf} -ERC GVCF\n"
    command << "gzip #{g_vcf}\n"
    command << "echo '#{GlobalVariables::SUSHI}' > #{@dataset['Name']}_dummy.txt\n"
  end
end

if __FILE__ == $0

end
