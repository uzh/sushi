#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class GatkJoinGenoTypesRNASeqApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'GATK JoinGenotypesRNASeq'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Variants'
    @description =<<-EOS
Joint genotyping of RNA-seq gVCF files (GATK GenotypeGVCFs), followed by<br/>
- hard filtering after the GATK RNA-seq best practices (VariantFiltration: FS > 30, QD < 2, clusters of 3 SNVs within 35 bp; FILTER column set to PASS or the filter name, plus a PASS-only VCF)<br/>
- splitting of multiallelic records and left-alignment of indels (bcftools norm)<br/>
- variant effect and impact annotation with SnpEff (database built from the reference annotation on first use)<br/>
    EOS
    @required_columns = ['Name','GVCF','GVCFINDEX','Species','refBuild']
    @required_params = ['name','grouping']
    @params['cores'] = '4'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '50'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '100'
    @params['scratch', "context"] = "slurm"
    @params['name'] = 'GATK_GenotypingRNASeq'
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "reference genome assembly"
    @params['grouping'] = ''
    @params['grouping', "context"] = "GatkJoinGenoTypesRNASeq"
    @params['minReadDepth'] = '20'
    @params['minReadDepth', 'description'] = 'genotypes with a lower read depth are shown as low coverage in the report'
    @params['hardFilter'] = true
    @params['hardFilter', 'description'] = 'GATK RNA-seq hard filters: FS > 30, QD < 2, SNP clusters (3 within 35 bp); writes FILTER=PASS or the filter name, plus a PASS-only VCF'
    @params['normalize'] = true
    @params['normalize', 'description'] = 'split multiallelic records and left-align indels with bcftools norm'
    @params['snpEff'] = true
    @params['snpEff', 'description'] = 'annotate variant effects and impact (ANN field) with SnpEff; the first run for a reference builds the database (~15-30 min)'
    @params['specialOptions'] = ''
    @params['specialOptions', 'description'] = 'additional options passed to GATK GenotypeGVCFs'
    @params['mail'] = ""
    @modules = ["Dev/jdk", "Variants/GATK", "Tools/bcftools", "Tools/Picard", "Variants/SnpEff", "Dev/R"]
    @inherit_columns = ["Order Id"]
  end
  def next_dataset
    report_dir = File.join(@result_dir, @params['name'])
    dataset = {'Name'=>@params['name'],
     'VCF [File]'=>File.join(report_dir, "#{@params['name']}.vcf.gz"),
     'TBI [File]'=>File.join(report_dir, "#{@params['name']}.vcf.gz.tbi")
    }
    if @params['hardFilter'].to_s == 'true'
      dataset['PASS VCF [File]'] = File.join(report_dir, "#{@params['name']}.PASS.vcf.gz")
      dataset['PASS TBI [File]'] = File.join(report_dir, "#{@params['name']}.PASS.vcf.gz.tbi")
    end
    dataset.merge({
     'Report [File]'=>report_dir,
     'Html [Link]'=>File.join(report_dir, '00index.html'),
     'Species'=>(first = @dataset.first and first['Species']),
     'refBuild'=>@params['refBuild']
    }).merge(extract_columns(colnames: @inherit_columns))
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
  end

  def commands
   run_RApp("EzAppJoinGenoTypesRNASeq")
  end
      end
