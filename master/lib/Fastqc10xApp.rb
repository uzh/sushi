#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class Fastqc10xApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'Fastqc10x'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'QC'
    @description =<<-EOS
A quality control tool for NGS reads<br/>
<a href='http://www.bioinformatics.babraham.ac.uk/projects/fastqc'/>Web-site with docu and a tutorial video</a>
EOS
    @required_columns = ['Name','RawDataDir']
    @required_params = ['name', 'paired']
    @params['cores'] = '8'
    @params['ram'] = '16'
    @params['scratch'] = '300'
    @params['paired'] = true
    @params['name'] = 'FastQC_Result'
    @params['cmdOptions'] = ""
    @params['mail'] = ""
    @modules = ["QC/FastQC", "Tools/Picard", "Tools/samtools", "Tools/sambamba", "Dev/Python"]
  end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, '00index.html')
    reportMultiQC_link = File.join(report_file, 'multiqc_report.html')
    {'Name'=>@params['name'],
     'Report [File]'=>report_file,
     'Html [Link]'=>report_link,
     'MultiQC [Link]'=>reportMultiQC_link,
    }
  end
  def commands
    run_RApp("EzAppFastqc_10x")
  end
end
