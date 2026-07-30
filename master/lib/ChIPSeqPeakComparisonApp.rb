#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class ChIPSeqPeakComparisonApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'ChIPSeqPeakComparison'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'GeneRegulation'
    @description =<<-EOS
Compares annotated ChIP-seq peaks across samples: pairwise overlap and consensus
peaks, ENCODE-style QC metrics (FRiP, NSC/RSC, NRF/PBC1/PBC2, BigWig correlation),
and the signal profiles and genomic footprint of the top N peaks per sample.
Expects already-called peaks (BED / narrowPeak / MACS xls / xlsx) plus BAM and/or
BigWig tracks; the peak, BAM and BigWig columns are auto-detected and each report
section degrades gracefully when its inputs are absent. Read-only with respect to
the input dataset.<br/>
EOS
    @required_columns = ['Name']
    @required_params  = ['refBuild']
    @params['cores'] = '8'
    @params['cores', 'context'] = 'slurm'
    @params['ram'] = '30'
    @params['ram', 'context'] = 'slurm'
    @params['scratch'] = '100'
    @params['scratch', 'context'] = 'slurm'
    @params['refBuild'] = ref_selector
    @params['refBuild', 'context'] = 'reference genome assembly'
    @params['paired'] = false
    @params['paired', 'description'] = 'are the BAM libraries paired-end?'
    @params['markType'] = ['narrow', 'broad']
    @params['markType', 'description'] = 'narrow (TF) or broad (histone); drives QC thresholds'
    @params['topN'] = 5
    @params['topN', 'description'] = 'number of most-significant peaks (large log2 over background) to show per-sample coverage tracks for'
    @params['rankBy'] = ['signalValue', 'qValue', 'score', 'width']
    @params['minSamplesForConsensus'] = 1
    @params['minOverlapBp'] = 1
    @params['profileExtend'] = 2000
    @params['profileBinSize'] = 50
    @params['quantifyFrom'] = ['bigwig', 'bam']
    @params['normalization'] = ['CPM', 'none', 'quantile']
    @params['peakFormat'] = ['auto', 'narrowPeak', 'broadPeak', 'bed', 'macsXls', 'xlsx']
    @params['useBlacklist'] = true
    @params['blacklistFile'] = ''
    @params['blacklistFile', 'description'] = 'blacklist BED path; empty = none (e.g. non-model organisms)'
    @params['runDifferentialBinding'] = false
    @params['runDifferentialBinding', 'description'] = 'stub, off; the consensus count matrix is the input for a follow-up app'
    @params['name'] = 'ChIPSeqPeakComparison'
    @params['specialOptions'] = ''
    @params['mail'] = ''
    @modules = ['Dev/R']
    @inherit_tags = ['Factor', 'B-Fabric', 'Characteristic']
  end
  def next_dataset
    report_file = File.join(@result_dir, "#{@params['name']}")
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@params['name'],
     'Species'=>(dataset = @dataset.first and dataset['Species']),
     'refBuild'=>@params['refBuild'],
     'Report [File]'=>report_file,
     'Static Report [Link]'=>report_link,
    }.merge(extract_columns(@inherit_tags))
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild'] if dataset_has_column?('refBuild')
    if dataset_has_column?('paired')
      @params['paired'] = @dataset[0]['paired']
    end
  end
  def commands
    run_RApp('EzAppChIPSeqPeakComparison')
  end
end

if __FILE__ == $0
end
