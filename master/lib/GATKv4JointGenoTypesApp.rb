#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class GATKv4JointGenoTypesApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'GATKv4JointGenoTypes'
    @params['process_mode'] = 'DATASET'
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
    @params['minDP'] = '2'
    @params['maxDP'] = '10000'
    @params['MQRankSum'] = '-15.0'
    @params['specialOptions'] = ''
    @modules = ["Variants/GATK/4.2.0.0", "Tools/Picard/2.18.0"]
    @inherit_columns = ["Order Id"]
  end
  def next_dataset
    {'Name'=>@params['name'],
     'Raw VCF [File]'=>File.join(@result_dir, "#{@params['name']}.raw.vcf.gz"),
     'Filtered VCF [File]'=>File.join(@result_dir, "#{@params['name']}.filtered.vcf.gz"),
     'Species'=>(dataset = @dataset.first and dataset['Species']),
     'refBuild'=>@params['refBuild']
    }.merge(extract_columns(colnames: @inherit_columns))
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
  end
  def commands
    refs = GENOME_REF_DIRS.map{|genome_ref_dir| File.join(genome_ref_dir, @params['refBuild'].split('/')[0,3].join("/")+"/Sequence/WholeGenomeFasta/genome.fa")}
    ref = refs.find{|fa| File.exist?(fa)}
    combined_g_vcf = @params['name'] + ".g.vcf"
    combined_raw_vcf = @params['name'] + ".raw.vcf"
    combined_filter_marked_vcf = @params['name'] + ".filter_marked.vcf"
    combined_filtered_vcf = @params['name'] + ".filtered.vcf"
    gvcfs = []
    command = ""
    @dataset.each do |row|
      gvcf_gz = File.join(@gstore_dir, row['GVCF [File]'])
      gvcf = File.basename(gvcf_gz).gsub(/.gz/, '')
      gvcfs << gvcf
      command << "cp #{gvcf_gz} .\n"
      command << "gunzip -c #{File.basename(gvcf_gz)} > #{gvcf}\n"
      command << "rm #{File.basename(gvcf_gz)}\n"
    end
    gvcfs_option = gvcfs.map{|gvcf| "-V #{gvcf}"}.join(" ")
    command << "gatk CombineGVCFs -R #{ref} #{gvcfs_option} -O #{combined_g_vcf}\n"
    command << "gatk GenotypeGVCFs -R #{ref} -V #{combined_g_vcf} -O #{combined_raw_vcf}\n"
    #raw_vcf = combined_raw_vcf
    #command << "gatk VariantFiltration -R #{ref} -V #{raw_vcf} --filter-expression \"! vc.hasAttribute('QD') || QD < #{@params['QD']}\" --filter-name \"QD\" --filter-expression \"vc.isSNP() && (MQ < #{@params['MQ']} || (vc.hasAttribute('MQRankSum') && MQRankSum < #{@params['MQRankSum']}))\" --filter-name \"MQ\" --genotype-filter-expression \"GQ < #{@params['GQ']} || DP == #{@params['DP']}\" --genotype-filter-name \"GQ\" -O #{combined_filtered_vcf}\n"
    command << "gatk VariantFiltration -R #{ref} -V #{combined_raw_vcf} --filter-name \"QD\" --filter-expression \"vc.hasAttribute('QD') && QD < #{@params['QD']}\" --filter-name \"MQ\" --filter-expression \"vc.isSNP() && vc.hasAttribute('MQ') && MQ < #{@params['MQ']}\" --filter-name \"MQRankSum\" --filter-expression \"vc.hasAttribute('MQRankSum') && MQRankSum < #{@params['MQRankSum']}\" --genotype-filter-name \"GQ\" --genotype-filter-expression \"GQ < #{@params['GQ']}\" --filter-name \"DP\" --filter-expression \"vc.hasAttribute('DP') && (DP < #{@params['minDP']} || DP >= #{@params['maxDP']})\" -O #{combined_filter_marked_vcf}\n"
    if @params['only_SNP']
      command << "gatk SelectVariants -R #{ref} -V #{combined_filter_marked_vcf} -O #{combined_filtered_vcf} -select-type SNP --exclude-filtered\n"
    else
      command << "gatk SelectVariants -R #{ref} -V #{combined_filter_marked_vcf} -O #{combined_filtered_vcf} --exclude-filtered\n"
    end
    command << "gzip #{combined_raw_vcf}\n"
    command << "gzip #{combined_filtered_vcf}\n"
    command
  end
end
