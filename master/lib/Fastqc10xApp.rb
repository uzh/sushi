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
    @required_columns = ['Name','RawDataDir', 'Read Count']
    @required_params = ['name', 'paired']
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '30'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '300'
    @params['scratch', "context"] = "slurm"
    @params['paired'] = true
    @params['paired', "context"] = "Fastqc10x"
    @params['name'] = 'FastQC_Result'
    @params['cmdOptions'] = ""
    @params['cmdOptions', "context"] = "Fastqc10x"
    @params['mail'] = ""
    @modules = ["QC/FastQC", "Tools/Picard", "Tools/samtools", "Dev/Python"]
    @inherit_columns = ["Order Id"]
  end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, '00index.html')
    reportMultiQC_link = File.join(report_file, 'multiqc_report.html')
    {'Name'=>@params['name'],
     'Report [File]'=>report_file,
     'Html [Link]'=>report_link,
     'MultiQC [Link]'=>reportMultiQC_link,
    }.merge(extract_columns(colnames: @inherit_columns))
  end
  def commands
    run_RApp("EzAppFastqc_10x")
  end
end
