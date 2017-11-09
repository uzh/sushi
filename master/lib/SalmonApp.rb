#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class SalmonApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Salmon'
    @analysis_category = 'Count'
    @description =<<-EOS
	 <a href="https://github.com/COMBINE-lab/salmon">Salmon</a> is a wicked-fast program to produce a highly-accurate, transcript-level quantification estimates from RNA-seq data.
EOS
    @required_columns = ['Name','Read1','Species']
    @required_params = ['paired', 'strandMode']
    @params['paired'] = false
    @params['strandMode'] = ['both', 'sense', 'antisense']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '100'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['fldMean'] = 0
    @params['fldMean', 'description'] = 'expected mean fragment length for single-end reads'
    @params['fldSD'] = 0
    @params['fldSD', 'description'] = 'expected standard deviation of the fragment length distribution for single-end reads'
    @params['useVBOpt'] = false
    @params['useVBOpt', 'description'] = 'use the variational Bayesian EM algorithm rather than the "standard" EM algorithm to optimize abundance estimates'
    @params['numBootstraps'] = 0
    @params['numBootstraps', 'description'] = 'number of bootstrapped samples; by default, Salmon does not use bootstrapping'
    @params['numGibbsSamples'] = 0
    @params['numGibbsSamples', 'description'] = 'an alternative to bootstrapping for estimating the variance in abundance estimates'
    @params['seqBias'] = false
    @params['seqBias', 'description'] = 'enable sequence bias correction'
    @params['gcBias'] = false
    @params['gcBias', 'description'] = 'enable GC bias correction'
    @params['posBias'] = false
    @params['posBias', 'description'] = 'enable position bias correction'
    @params['specialParams'] = ''
    @params['specialParams', 'description'] = 'additional command line parameters to pass to Salmon'
    @params['transcriptFasta'] = ''
    @params['transcriptFasta', 'description'] = 'give full path of transcript fasta file; in that case the build is ignored; if it comes from trinity assembly the gene-isoform associations will be extracted and used'
    @params['transcriptTypes'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA', 'long_noncoding', 'short_noncoding', 'pseudogene']
    @params['transcriptTypes', 'multi_selection'] = true
    @params['trimAdapter'] = true
    @params['trimLeft'] = 0
    @params['trimRight'] = 0
    @params['minTailQuality'] = 0
    @params['minAvgQuality'] = 20
    @params['mail'] = ""
    @modules = ["Aligner/Salmon", "QC/Flexbar", "QC/Trimmomatic", "Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    else
      @required_columns << 'fldMean' << 'fldSD'
    end
  end
  def set_default_parameters
    @params['paired'] = dataset_has_column?('Read2')
    if dataset_has_column?('strandMode')
      @params['strandMode'] = @dataset[0]['strandMode']
    end
  end
  def next_dataset
    {
      'Name'=>@dataset['Name'],
      'Count [File]'=>File.join(@result_dir, "#{@dataset['Name']}.txt"),
      'Species'=>@dataset['Species'],
      'refBuild'=>@params['refBuild'],
      'featureLevel'=>'isoform',
      'refFeatureFile'=>@params['refFeatureFile'],
      'strandMode'=>@params['strandMode'],
      'paired'=>@params['paired'],
      'Read Count'=>@dataset['Read Count'],
      'transcriptTypes'=>@params['transcriptTypes']
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppSalmon")
  end
end

if __FILE__ == $0
  run SalmonApp

end
