#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class SCRNAVelocityApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'SCRNAVelocity'
    @params['process_mode'] = 'SAMPLE'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
RNA velocity analysis for single cell data<br/>
    EOS
    @required_columns = ['Name', 'refBuild', 'Report', 'ResultDir']
    @required_params = ['name']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '16'
    @params['scratch'] = '150'
    @params['name'] = 'SCRNAVelocity'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['scProtocol'] = ['smart-Seq2', '10X']
    @params['markersToCheck'] = ''
    @params['markersToCheck', 'description'] = 'Show the fitting of individual genes if given, in format of Lgals3,Napsa.'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Dev/R", "Dev/Python/3.6.4", "Tools/samtools"]
  end
  def next_dataset
    report_file = File.join(@result_dir, "#{@dataset['Name']}_SCRNAVelocity")
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@dataset['Name'],
     'refBuild'=>@params['refBuild'],
     'refFeatureFile'=>@params['refFeatureFile'],
     'Static Report [Link]'=>report_link,
     'Report [File]'=>report_file,
    }
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('refFeatureFile')
      @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    end
  end
  def commands
    run_RApp("EzAppSCRNAVelocity")
  end
end

if __FILE__ == $0

end

