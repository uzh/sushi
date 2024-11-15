#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class GATKJointGenotypingByGenomicsDBApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'GATKJointGenotypingByGenomicsDBApp'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Variants'
    @description =<<-EOS
Genotyping by GenomicsDBImport,GenotypeGVCFs, and hard-filtering by VariantFiltration 
genotype,merge and annotate gvcf-Files<br/>
    EOS
    @required_columns = ['Name', 'GVCF', 'GVCFINDEX', 'Species', 'refBuild', 'Dummy']
    @required_params = ['name']
    @params['cores'] = '8'
    @params['ram'] = '50'
    @params['scratch'] = '100'
    @params['name'] = 'GATKJointGenotyping'
    @params['refBuild'] = ref_selector

    @params['ConcGCThreads', 'hr-header'] = "JavaVM options"
    @params['ConcGCThreads'] = '1'
    @params['ConcGCThreads', 'description'] = 'Number of threads for releasing memory'
    @params['ParallelGCThreads'] = '8'
    @params['ParallelGCThreads', 'description'] = 'Number of threads concurrent garbage collectors will use'
    @params['Xmx'] = '50'
    @params['Xmx', 'description'] = 'Maximum heap memory size for JVM'

    #@params['batch-size', 'hr-header'] = 'GenomicsDBImport options'
    #@params['batch-size'] = '10'
    #@params['reader-threads'] = '1'
    @params['QD', 'hr-header'] = 'VariantFiltration options'
    @params['QD'] = '2.0'
    @params['FS'] = '60.0'
    @params['MQ'] = '40.0'
    @params['GQ'] = '20'
    #@params['ReadPosRankSum'] = '-8.0'
    @params['MQRankSum'] = '-12.5'
    @modules = ["Variants/GATK/4.5.0.0"]
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
    annos = GENOME_REF_DIRS.map{|genome_ref_dir| File.join(genome_ref_dir, @params['refBuild']+"/Genes/intervals.list")}
    intervals_list = annos.find{|list| File.exist?(list)}
    combined_raw_vcf = @params['name'] + ".raw.vcf.gz"
    combined_filtered_vcf = @params['name'] + ".filtered.vcf.gz"
		jvm_options = "-XX:ConcGCThreads=#{@params["ConcGCThreads"]} -Xmx#{@params['Xmx']}G -XX:ParallelGCThreads=#{@params['ParallelGCThreads']}"
		command = <<EOS
> sample_list.txt
tail -n +2 "\$INPUT_DATASET" | while IFS=$'\\t' read -r name gvcf gvcfindex rest; do
    gvcf_path="\$GSTORE_DIR/\$gvcf"
    echo -e "\$name\\t\$gvcf_path" >> sample_list.txt
done
EOS
    command << if intervals_list
								 "gatk --java-options #{jvm_options} GenomicsDBImport -R #{ref} --sample-name-map sample_list.txt -L #{intervals_list} --genomicsdb-workspace-path genomics_db\n"
							 else
								 "gatk --java-options #{jvm_options} GenomicsDBImport -R #{ref} --sample-name-map sample_list.txt --genomicsdb-workspace-path genomics_db\n"
							 end
    command << "gatk --java-options #{jvm_options} GenotypeGVCFs -R #{ref} -V gendb://genomics_db -O #{combined_raw_vcf}\n"
    command +=<<EOS
gatk --java-options #{jvm_options} VariantFiltration -R #{ref} -V #{combined_raw_vcf} -O #{combined_filtered_vcf} \\
--filter-name "QD" --filter-expression "vc.hasAttribute('QD') && QD < #{@params['QD']}" \\
--filter-name "FS" --filter-expression "FS > #{@params['FS']}" \\
--filter-name "MQ" --filter-expression "vc.isSNP() && vc.hasAttribute('MQ') && MQ < #{@params['MQ']}" \\
--filter-name "MQRankSum" --filter-expression "vc.hasAttribute('MQRankSum') && MQRankSum < #{@params['MQRankSum']}" \\
--genotype-filter-name "GQ" --genotype-filter-expression "GQ < #{@params['GQ']}"
EOS
    command
  end
end
