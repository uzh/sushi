#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class NcPROApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'ncPRO_Report'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Count'
    @description =<<-EOS
Annotation and Profiling of ncRNAs in smallRNA-seq<br/>
Uses ncPRO-seq for a complete analysis of small-RNA-seq. ncPRO-seq considers performs quality assessment and quantiation.
It considers various species of small RNA.<br />
<a href='https://ncpro.curie.fr/index.html'>https://ncpro.curie.fr/index.html/</a>
ncPRO can only process single-end stranded RNA.
<strong>IMPORTANT: ncPRO does not run on the virtual nodes!!!!</strong>
EOS

    @required_columns = ['Name','Read1', 'Adapter1', 'Species']
    @required_params = ['name', 'cores', 'ram', 'scratch']

    @params['cores'] = '8'
    @params['ram'] = '40'
    @params['scratch'] = '100'
    @params['refBuild'] = ['Mus_musculus/UCSC/mm10', 'Homo_sapiens/UCSC/hg19', 'Rattus_norvegicus/UCSC/rn5']
    #@params['paired'] = false
    #@params['paired', 'description'] = 'whether the reads are paired end; must be false since ncPRO-seq does not support paired-end'
    @params['name'] = 'ncPRO_Result'
    @params['mail'] = ""
    @modules = ["Aligner/Bowtie", "QC/Flexbar", "Tools/ncPROseq", "QC/Trimmomatic", "Dev/R"]
  end
  def preprocess
    if @params['paired']
      @required_columns<<  'Read2'
    end
  end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, 'ncpro/report.html')
    {'Name'=>@params['name'],
     'Species'=>(dataset = @dataset.first and dataset['Species']),
     'refBuild'=>@params['refBuild'],
     'Report [File]'=>report_file,
     'Html [Link]'=>report_link,
     'TrimCounts [Link]'=>File.join(report_file, 'trimCounts-barplot.png')
    }
  end
  def commands
    run_RApp("EzAppNcpro")
  end
end

if __FILE__ == $0
  run NcPROApp

end
