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
    @required_columns = ['Name','Read1','Species','PlateName']
    @required_params = ['name']

@params['cores'] = '8'
@params['ram'] = '20'
@params['scratch'] = '100'
@params['name'] = 'RNA_Bias_Result'
@params['paired'] = false
@params['strandMode'] = ['both', 'sense', 'antisense']
@params['refBuild'] = ref_selector
@params['refFeatureFile'] = 'genes.gtf'
@params['bootstrap-samples'] = '10'
@params['seed'] = '42'
@params['fragment-length'] = '150'
@params['sd'] = '70'
@params['bias'] = true
@params['pseudobam'] = true
@params['transcriptFasta'] = ''
@params['transcriptTypes'] = 'protein_coding'
@params['trimAdapter'] = true
@params['cut_front'] = false
@params['trim_front1'] = '0'
@params['trim_tail1'] = '0'
@params['cut_front'] = false
@params['cut_front_window_size'] =	'4'
@params['cut_front_mean_quality'] =	'20'
@params['cut_tail']	= true
@params['cut_tail_window_size'] =	'4'
@params['cut_tail_mean_quality'] =	'15'
@params['cut_right'] =	false
@params['cut_right_window_size'] =	'4'
@params['cut_right_mean_quality'] =	'20'
@params['average_qual'] =	'20'
@params['max_len1'] = 	'0'
@params['max_len2'] =	'0'
@params['poly_x_min_len'] =	'10'
@params['length_required'] =	'20'
@params['backgroundExpression'] = 5 ## low numbers needed for the iSeq runs with
@params['sigThresh'] = 5
@params['normMethod'] = 'logMean'
@params['expressionName'] = 'est_counts'
@params['minReadsPerSample']= '10000'
@params['cmdOptions'] = ''
@params['mail'] = ''

@modules = ["Dev/R", "Aligner/kallisto/0.46.1_deb10", "QC/fastp", "Tools/samtools"]
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
     'Report [File]'=>report_file,
     'Report [Link]'=>report_file
    }.merge(extract_columns(colnames: @inherit_columns))
  end
  def commands
    run_RApp("EzAppRnaComputeBias")
  end
end

if __FILE__ == $0
  
end
