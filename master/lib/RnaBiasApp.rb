#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class RnaBiasApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'RnaBiasApp'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'QC'
    @description =<<-EOS 
EOS
    @required_columns = ['Name','Read1','Species','Order Id','PlateName']
    @required_params = ['name']

@params['cores'] = '8'
@params['cores', "context"] = "slurm"
@params['ram'] = '50'
@params['ram', "context"] = "slurm"
@params['scratch'] = '100'
@params['scratch', "context"] = "slurm"
@params['name'] = 'RNA_Bias_Result'
@params['paired'] = false
@params['paired', "context"] = "RnaBias"
@params['strandMode'] = ['both', 'sense', 'antisense']
@params['strandMode', "context"] = "RnaBias"
@params['refBuild'] = ref_selector
@params['refBuild', "context"] = "referfence genome assembly"
@params['refFeatureFile'] = 'genes.gtf'
@params['refFeatureFile', "context"] = "RnaBias"
@params['bootstrap-samples'] = '10'
@params['seed'] = '42'
@params['fragment-length'] = '150'
@params['sd'] = '70'
@params['transcriptFasta'] = ''
@params['transcriptTypes'] = 'protein_coding'
@params['trimAdapter'] = true
@params['trimAdapter', "context"] = "OpenGene/fastp"
@params['cut_front'] = false
@params['cut_front', "context"] = "OpenGene/fastp"
@params['trim_front1'] = '0'
@params['trim_front1', "context"] = "OpenGene/fastp"
@params['trim_tail1'] = '0'
@params['trim_tail1', "context"] = "OpenGene/fastp"
@params['cut_front'] = false
@params['cut_front_window_size'] =	'4'
@params['cut_front_window_size', "context"] = "OpenGene/fastp"
@params['cut_front_mean_quality'] =	'20'
@params['cut_front_mean_quality', "context"] = "OpenGene/fastp"
@params['cut_tail']	= true
@params['cut_tail', "context"] = "OpenGene/fastp"
@params['cut_tail_window_size'] =	'4'
@params['cut_tail_window_size', "context"] = "OpenGene/fastp"
@params['cut_tail_mean_quality'] =	'15'
@params['cut_tail_mean_quality', "context"] = "OpenGene/fastp"
@params['cut_right'] =	false
@params['cut_right', "context"] = "OpenGene/fastp"
@params['cut_right_window_size'] =	'4'
@params['cut_right_window_size', "context"] = "OpenGene/fastp"
@params['cut_right_mean_quality'] =	'20'
@params['cut_right_mean_quality', "context"] = "OpenGene/fastp"
@params['average_qual'] =	'20'
@params['average_qual', "context"] = "OpenGene/fastp"
@params['max_len1'] = 	'0'
@params['max_len1', "context"] = "OpenGene/fastp"
@params['max_len2'] =	'0'
@params['max_len2', "context"] = "OpenGene/fastp"
@params['poly_x_min_len'] =	'10'
@params['poly_x_min_len', "context"] = "OpenGene/fastp"
@params['length_required'] =	'20'
@params['length_required', "context"] = "OpenGene/fastp"
@params['backgroundExpression'] = 5 ## low numbers needed for the iSeq runs with
@params['backgroundExpression', "context"] = "RnaBias"
@params['sigThresh'] = 5
@params['sigThresh', "context"] = "RnaBias"
@params['normMethod'] = 'logMean'
@params['normMethod', "context"] = "RnaBias"
@params['expressionName'] = 'est_counts'
@params['expressionName', "context"] = "RnaBias"
@params['minReadsPerSample']= '10000'
@params['minReadsPerSample', "context"] = "RnaBias"
@params['cmdOptions'] = ''
@params['cmdOptions', "context"] = "RnaBias"
@params['mail'] = ''

@modules = ["Dev/R", "Aligner/kallisto", "QC/fastp", "Tools/samtools", "Aligner/STAR", "Dev/Python", "Tools/Picard", "Dev/jdk"]
@inherit_columns = ["Order Id"]
  end
 def set_default_parameters
    @params['paired'] = dataset_has_column?('Read2')
    if dataset_has_column?('strandMode')
      @params['strandMode'] = @dataset[0]['strandMode']
    end
 end
  def preprocess
    if @params['paired']
      @required_columns<<  'Read2'
    end
  end
  def next_dataset
    report_file = File.join(@result_dir,'report-RNAseq-bias.pdf')    
    {'Name'=>@params['name'],
     'Report [File,Link]'=>report_file,
    }.merge(extract_columns(colnames: @inherit_columns))
  end
  def commands
    run_RApp("EzAppRnaComputeBias")
  end
end

if __FILE__ == $0
  
end
