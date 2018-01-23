#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-095657'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class TophatApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Tophat'
    @analysis_category = 'Map'
    @description =<<-EOS
A spliced read mapper for RNA-Seq<br/>
<a href='https://ccb.jhu.edu/software/tophat/index.shtml'>manual</a>
    EOS
    @required_columns = ['Name','Read1','Species']
    @required_params = ['refBuild','paired', 'strandMode']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '16'
    @params['scratch'] = '100'
    @params['refBuild'] = ref_selector
    @params['paired'] = false
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['refFeatureFile'] = 'genes.gtf'
    @params['cmdOptions'] = '--mate-inner-dist 100 --mate-std-dev 150'
    @params['trimAdapter'] = true
    @params['trimLeft'] = 0
    @params['trimRight'] = 0
    @params['minTailQuality'] = 0
    @params['minAvgQuality'] = 10
    @params['minReadLength'] = 20
    @params['specialOptions'] = ''
    @params['mail'] = ""
    # Python2 is required because of RSeQC package
    @modules = ["Tools/samtools", "Aligner/Bowtie", "Aligner/Bowtie2", "Aligner/TopHat", "QC/Flexbar", "QC/Trimmomatic", "Dev/Python2", "Dev/R", "Tools/sambamba"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  end
  def set_default_parameters
    @params['paired'] = dataset_has_column?('Read2')
    if dataset_has_column?('strandMode')
      @params['strandMode'] = @dataset[0]['strandMode']
    end
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'BAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam"),
     'BAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam.bai"),
     'IGV Starter [Link]'=>File.join(@result_dir, "#{@dataset['Name']}-igv.jnlp"),
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'paired'=>@params['paired'],
     'refFeatureFile'=>@params['refFeatureFile'],
     'strandMode'=>@params['strandMode'],
     'Read Count'=>@dataset['Read Count'],
     'IGV Starter [File]'=>File.join(@result_dir, "#{@dataset['Name']}-igv.jnlp"),
     'IGV Session [File]'=>File.join(@result_dir, "#{@dataset['Name']}-igv.xml"),
     'PreprocessingLog [File]'=>File.join(@result_dir, "#{@dataset['Name']}_preprocessing.log")
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppTophat")
  end
end

if __FILE__ == $0
  run TophatApp
end
