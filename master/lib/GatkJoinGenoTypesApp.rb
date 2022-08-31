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
    @params['specialOptions'] = ''
    @params['dbNSFP_file'] = ["dbNSFP2.9.txt.gz", "dbNSFP4.1c.txt.gz"]
    @params['dbNSFP_fields'] = 'ExAC_AF,ExAC_AC,1000Gp1_EUR_AF,Uniprot_acc,Interpro_domain,phastCons100way_vertebrate,CADD_phred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,SIFT_score,SIFT_pred'
    @params['mail'] = ""
    @modules = ["Dev/jdk", "Variants/GATK", "Tools/Picard", "Variants/SnpEff", "Dev/R"]
  end
  def next_dataset
    report_dir = File.join(@result_dir, @params['name'])
    {'Name'=>@params['name'],
     'Report [File]'=>report_dir,
#     'Html [Link]'=>File.join(report_dir, '00index.html'),
     'Species'=>(dataset = @dataset.first and dataset['Species']),
     'refBuild'=>@params['refBuild']
    }
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
