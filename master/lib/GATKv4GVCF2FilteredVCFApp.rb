#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class GATKv4GVCF2FilteredVCFApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'GATKv4GVCF2FilteredVCF'
    @analysis_category = 'Variants'
    @description =<<-EOS
genotype,merge and annotate gvcf-Files<br/>
    EOS
    @required_columns = ['Name', 'GVCF', 'GVCFINDEX', 'Species', 'refBuild', 'Dummy']
    @required_params = ['name']
    @params['cores'] = '1'
    @params['ram'] = '50'
    @params['scratch'] = '100'
    @params['name'] = 'GATKv4_Genotyping'
    @params['refBuild'] = ref_selector
    @params['only_SNP'] = true
    @params['QD'] = '2.0'
    @params['MQ'] = '30.0'
    @params['GQ'] = '20'
    @params['DP'] = '0'
    @params['MQRankSum'] = '-15.0'
    @params['specialOptions'] = ''
    @modules = ["Variants/GATK/4.1.2.0", "Tools/Picard/2.18.0"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'Raw VCF [File]'=>File.join(@result_dir, "#{@dataset['Name']}.raw.vcf.gz"),
     'Filtered VCF [File]'=>File.join(@result_dir, "#{@dataset['Name']}.filtered.vcf.gz"),
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild']
    }
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
  end
  def commands
    ref = File.join(GENOME_REF_DIR, @params['refBuild'].split('/')[0,3].join("/")+"/Sequence/WholeGenomeFasta/genome.fa")
    combined_raw_vcf = @dataset['Name'] + ".raw.vcf"
    combined_raw_snp_vcf = @dataset['Name'] + ".raw.snp.vcf"
    combined_filtered_vcf = @dataset['Name'] + ".filtered.vcf"
    command = ""
    gvcf_gz = File.join(@gstore_dir, @dataset['GVCF'])
    gvcf = File.basename(gvcf_gz).gsub(/.gz/, '')
    command << "cp #{gvcf_gz} .\n"
    command << "gunzip -c #{File.basename(gvcf_gz)} > #{gvcf}\n"
    command << "rm #{File.basename(gvcf_gz)}\n"

    command << "gatk GenotypeGVCFs -R #{ref} -V #{gvcf} -O #{combined_raw_vcf}\n"
    raw_vcf = combined_raw_vcf
    if @params['only_SNP']
      command << "gatk SelectVariants -R #{ref} -V #{combined_raw_vcf} -O #{combined_raw_snp_vcf} -select-type SNP\n"
      raw_vcf = combined_raw_snp_vcf
    end
    command << "gatk VariantFiltration -R #{ref} -V #{raw_vcf} --filter-expression \"! vc.hasAttribute('QD') || QD < #{@params['QD']}\" --filter-name \"QD\" --filter-expression \"vc.isSNP() && (MQ < #{@params['MQ']} || (vc.hasAttribute('MQRankSum') && MQRankSum < #{@params['MQRankSum']}))\" --filter-name \"MQ\" --genotype-filter-expression \"GQ < #{@params['GQ']} || DP == #{@params['DP']}\" --genotype-filter-name \"GQ\" -O #{combined_filtered_vcf}\n"
    if @params['only_SNP']
      command << "mv #{combined_raw_snp_vcf} #{combined_raw_vcf}\n"
    end
    command << "gzip #{combined_raw_vcf}\n"
    command << "gzip #{combined_filtered_vcf}\n"
    command
  end
end
