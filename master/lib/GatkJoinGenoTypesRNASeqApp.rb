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
genotype,merge and annotate gvcf-Files<br/>
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
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Dev/jdk", "Variants/GATK", "Tools/bcftools", "Tools/Picard", "Dev/R"]
    @inherit_columns = ["Order Id"]
  end
  def next_dataset
    report_dir = File.join(@result_dir, @params['name'])
    {'Name'=>@params['name'],
     'Report [File]'=>report_dir,
     'Html [Link]'=>File.join(report_dir, '00index.html'),
     'Species'=>(dataset = @dataset.first and dataset['Species']),
     'refBuild'=>@params['refBuild']
    }.merge(extract_columns(colnames: @inherit_columns))
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
  end

  def commands
   run_RApp("EzAppJoinGenoTypesRNASeq")
  end
      end
