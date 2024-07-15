#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class GatkJoinGenoTypesApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'GATK JoinGenotypes'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Variants'
    @description =<<-EOS
genotype,merge and annotate gvcf-Files<br/>
    EOS
    @required_columns = ['Name','GVCF','GVCFINDEX','Species','refBuild']
    @required_params = ['name','grouping']
    @params['cores'] = '8'
    @params['ram'] = '50'
    @params['scratch'] = '100'
    @params['name'] = 'GATK_Genotyping'
    @params['refBuild'] = ref_selector
    @params['grouping'] = ''
    @params['targetFile'] = ''
    @params['recalibrateVariants']=true
    @params['recalibrateVariants', 'description'] = "perform Variant Quality Score Recalibration (VQSR)"
    @params['recalibrateInDels']=false
    @params['recalibrateInDels', 'description'] = "perform InDel recalibration. This tends to crash for smaller or targeted datasets"
    @params['snpEffDB']=''
    @params['snpEffDB', 'description'] = "perform functional annotation with SNPEff if annotation database is installed"
    @params['proteinCodingTranscriptsOnly']=false
    @params['proteinCodingTranscriptsOnly', 'description'] = "SNPEff annotate variants only with protCod variants"
    @params['specialOptions'] = ''
    @params['dbNSFP_file'] = ["dbNSFP4.3_GRCh38.txt.gz"]
    @params['dbNSFP_file','description'] = 'human data specific database to add MAF and prediction scores. Please use the version which matches your reference.'
    @params['dbNSFP_fields'] = 'CADD_phred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,DANN_score,phastCons100way_vertebrate,gnomAD_exomes_AF,gnomAD_exomes_AC,gnomAD_genomes_AF,gnomAD_exomes_AC,clinvar_id,clinvar_clnsig,clinvar_review,clinvar_OMIM_id'
    @params['dbNSFP_fields','description'] = "details about available fields are available under <a href='https://sites.google.com/site/jpopgen/dbNSFP',>dbNSFP</a>"
    @params['mail'] = ""
    @modules = ["Dev/jdk", "Variants/GATK", "Tools/Picard", "Variants/SnpEff", "Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def next_dataset
    report_dir = File.join(@result_dir, @params['name'])
    {'Name'=>@params['name'],
     'Report [File]'=>report_dir,
#     'Html [Link]'=>File.join(report_dir, '00index.html'),
     'Species'=>(dataset = @dataset.first and dataset['Species']),
     'refBuild'=>@params['refBuild']
    }.merge(extract_columns(@inherit_tags))
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('targetFile')
      @params['targetFile'] = @dataset[0]['targetFile']
    end

  end

  def commands
   run_RApp("EzAppJoinGenoTypes")
  end
      end
