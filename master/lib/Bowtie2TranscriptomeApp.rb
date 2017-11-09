#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-093702'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class Bowtie2TranscriptomeApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Bowtie2TranscriptomeApp'
    @analysis_category = 'Map'
    @description =<<-EOS
Aligns the reads to a transcriptome only index. The index is built from the gtf file specified. For alignments to
genome or joined see the other mappers.
Mapping to transcriptome is useful, e.g. for Ribo-Seq data.
For expression estimation of RNA-seq data we recommend using the RSEM counting App or STAR+featureCounts.
For bowtie2 options see:<br/>
<a href='http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml'>http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml/</a>
EOS

    @required_columns = ['Name','Read1','Species']
    @required_params = ['refBuild','paired']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '16'
    @params['scratch'] = '100'
    @params['refBuild'] = ref_selector
    @params['refBuild', 'description'] = 'the genome refBuild and annotation to use as reference. If human variant calling is the main goal, please use hg_19_karyotypic.'
    @params['paired'] = false
    @params['paired', 'description'] = 'whether the reads are paired end; if false then only Read1 is considered even if Read2 is available.'
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['refFeatureFile'] = 'genes.gtf'
    @params['cmdOptions'] = '--no-unal'
    @params['cmdOptions', 'description'] = 'specify the commandline options for bowtie2; do not specify any option that is already covered by the dedicated input fields'
    @params['trimAdapter'] = true
    @params['trimAdapter', 'description'] = 'if adapters should be trimmed'
    @params['trimLeft'] = 0
    @params['trimLeft', 'description'] = 'fixed trimming at the "left" i.e. 5-prime end of the read'
    @params['trimRight'] = 0
    @params['trimRight', 'description'] = 'fixed trimming at the "right" i.e. 3-prime end of the read'
    @params['minTailQuality'] = 0
    @params['minTailQuality', 'description'] = 'if above zero, then reads are trimmed as soon as 4 consecutive bases have lower mean quality'
    @params['minAvgQuality'] = 10
    @params['minReadLength'] = 20
    @params['specialOptions'] = ''
    @params['specialOptions', 'description'] = 'special unsupported options that the R wrapper may support, format: <key>=<value>'
    @params['mail'] = ""
    @modules = ["Tools/samtools", "Aligner/Bowtie2", "QC/Flexbar", "QC/Trimmomatic", "Dev/R"]
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
     'trBAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam"),
     'trBAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam.bai"),
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'paired'=>@params['paired'],
     'refFeatureFile'=>@params['refFeatureFile'],
     'strandMode'=>@params['strandMode'],
     'Read Count'=>@dataset['Read Count']
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppBowtie2Transcriptome")
  end
end
