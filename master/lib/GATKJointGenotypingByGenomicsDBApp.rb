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
    @params['Xmx'] = '30'
    @params['Xmx', 'description'] = 'Maximum heap memory size for JVM'

    @params['intervals', 'hr-header'] = 'GenomicsDBImport options'
    @params['intervals'] = ''
    @params['intervals', 'file_upload'] = true
    @params['batch-size'] = '10'
    @params['reader-threads'] = '1'
    @params['QD', 'hr-header'] = 'VariantFiltration options'
    @params['QD'] = '2.0'
    @params['FS'] = '60.0'
    @params['MQ'] = '40.0'
    @params['GQ'] = '20'
    @params['ReadPosRankSum'] = '-8.0'
    @params['MQRankSum'] = '-12.5'
    @modules = ["Variants/GATK/4.6.0.0", "Tools/Picard/2.22.8"]
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
    combined_raw_vcf = @params['name'] + ".raw.vcf.gz"
    combined_filtered_vcf = @params['name'] + ".filtered.vcf.gz"
		jvm_options = ""
    command = ""
    sample_list = "sample_list.txt"
		intervals = "intervals.list"
		genomics_db = "genomics_db"
    open(sample_list, "w") do |f|
	    @dataset.each do |row|
	      name = row['Name']
	      gvcf_gz = File.join(@gstore_dir, row['GVCF [File]'])
				f.puts [name, gvcf_gz].join("\t")
	    end
    end
    command << "gatk #{jvm_options} GenomicsDBImport -R #{ref} --sample-name-map #{sample_list} -L #{intervals_list} --genomicsdb-workspace-path ${genomics_db}\n"
    command << "gatk GenotypeGVCFs -R #{ref} -V gendb://${genomics_db} -O #{combined_raw_vcf}\n"
    command << "gatk VariantFiltration -R #{ref} -V #{combined_raw_vcf} --filter-name \"QD\" --filter-expression \"vc.hasAttribute('QD') && QD < #{@params['QD']}\" --filter-name \"MQ\" --filter-expression \"vc.isSNP() && vc.hasAttribute('MQ') && MQ < #{@params['MQ']}\" --filter-name \"MQRankSum\" --filter-expression \"vc.hasAttribute('MQRankSum') && MQRankSum < #{@params['MQRankSum']}\" --genotype-filter-name \"GQ\" --genotype-filter-expression \"GQ < #{@params['GQ']}\" --filter-name \"DP\" --filter-expression \"vc.hasAttribute('DP') && (DP < #{@params['minDP']} || DP >= #{@params['maxDP']})\" -O #{combined_filter_vcf}\n"
    command
  end
end
