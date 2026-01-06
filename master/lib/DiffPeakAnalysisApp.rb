#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class DiffPeakAnalysisApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'DiffPeakAnalysis'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'GeneRegulation'
    @description =<<-EOS
    Finding differential peaks for ATAC Seq or ChIP Seq data <br/>
EOS
    @required_columns = ['Name','Count', 'BigWig']
    @required_params = ['grouping', 'sampleGroup', 'refGroup', 'refBuild']
    @params['cores'] = '4'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '20'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '100'
    @params['scratch', "context"] = "slurm"
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "reference genome assembly"
    @params['refFeatureFile'] = 'genes.gtf'
    @params['refFeatureFile', "context"] = "DiffPeakAnalysis"
    @params['grouping'] = ''
    @params['sampleGroup'] = ''
    @params['sampleGroup', 'description'] = 'sampleGroup should be different from refGroup'
    @params['refGroup'] = ''
    @params['refGroup', 'description'] = 'refGroup should be different from sampleGroup'
    @params['annotationMethod'] = ['chippeakanno', 'chipseeker', 'homer']
    @params['annotationMethod', 'description'] = 'peaks can be annotated with three different tools'
    @params['cmdOptions'] = ''
    @params['mail'] = ""
    @modules = ["Dev/R", "Tools/HOMER"]
    @inherit_columns = ["Order Id"]
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
  end
  
  def next_dataset
    @comparison = "#{@params['sampleGroup']}--over--#{@params['refGroup']}"
    @params['comparison'] = @comparison
    @params['name'] = @comparison
    report_file = File.join(@result_dir, "#{@params['comparison']}")
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@params['name'],
      'Report [Link]'=>report_link,
      'ResultFolder [File]'=>report_file
    }.merge(extract_columns(colnames: @inherit_columns))
  end
  def commands
    run_RApp("EzAppDiffPeakAnalysis")
  end
end

if __FILE__ == $0
  
end
