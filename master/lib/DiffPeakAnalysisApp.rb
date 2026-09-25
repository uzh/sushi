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
    Finding differential peaks for ATAC-Seq or ChIP-Seq data with DESeq2.<br/>
    Comparisons without biological replicates are supported (dispersion is
    estimated with a blind design and reported as conservative). The Quarto
    report includes volcano/MA/PCA/heatmap plots, peak feature and TSS-distance
    distributions, known-motif enrichment of up/down candidates (HOMER, against
    the consensus-peak background), GO over-representation with the peak-gene
    universe (mouse and human), BigWig signal profiles around the candidates,
    normalisation and GC/width bias diagnostics, one-click Enrichr submission of
    candidate genes, a downloadable result table, and session info.<br/>
EOS
    @required_columns = ['Name','Count']
    @required_params = ['grouping', 'sampleGroup', 'refGroup', 'refBuild']
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '30'
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
    @params['grouping2'] = ''
    @params['grouping2', 'description'] = 'optional second factor for paired/blocked testing (e.g. batch or donor). Type the dataset column name; the column must be tagged "NAME [Factor]". Leave empty for a single-factor analysis.'
    @params['grouping2', 'context'] = "DiffPeakAnalysis"
    @params['annotationMethod'] = ['homer', 'chippeakanno', 'chipseeker']
    @params['annotationMethod', 'description'] = 'peaks can be annotated with three different tools'
    @params['normMethod'] = ['DESeq2', 'readsInPeaks', 'TMM']
    @params['normMethod', 'description'] = 'size factors: DESeq2 median of ratios over peaks (assumes most peaks do not change), readsInPeaks (library size; use when many peaks change in one direction) or TMM. The report shows how the candidates would change under each method.'
    @params['normMethod', 'context'] = "DiffPeakAnalysis"
    @params['lfcThreshold'] = '1'
    @params['lfcThreshold', 'description'] = 'absolute log2 fold-change threshold for candidate peaks'
    @params['lfcThreshold', 'context'] = "DiffPeakAnalysis"
    @params['fdrThreshold'] = '0.05'
    @params['fdrThreshold', 'description'] = 'adjusted p-value threshold for candidate peaks'
    @params['fdrThreshold', 'context'] = "DiffPeakAnalysis"
    @params['lfcTest'] = false
    @params['lfcTest', 'description'] = 'test against |log2FC| > lfcThreshold instead of 0 (stricter, statistically proper fold-change cut)'
    @params['lfcTest', 'context'] = "DiffPeakAnalysis"
    @params['lfcShrink'] = ['apeglm', 'ashr', 'none']
    @params['lfcShrink', 'description'] = 'shrunken log2 fold change reported in an extra column (log2FoldChange_shrunk); candidates still use the unshrunken estimate'
    @params['lfcShrink', 'context'] = "DiffPeakAnalysis"
    @params['fitAllSamples'] = false
    @params['fitAllSamples', 'description'] = 'fit DESeq2 on all samples of the dataset (better dispersion estimates with many groups) and extract the sampleGroup vs refGroup contrast'
    @params['fitAllSamples', 'context'] = "DiffPeakAnalysis"
    @params['runMotifs'] = true
    @params['runMotifs', 'description'] = 'known-motif enrichment of up/down candidates with HOMER (any genome with a FASTA)'
    @params['runMotifs', 'context'] = "DiffPeakAnalysis"
    @params['motifDeNovo'] = false
    @params['motifDeNovo', 'description'] = 'also run HOMER de novo motif discovery (adds 30-60 min)'
    @params['motifDeNovo', 'context'] = "DiffPeakAnalysis"
    @params['runGoOra'] = true
    @params['runGoOra', 'description'] = 'GO biological-process over-representation of candidate genes (mouse and human only)'
    @params['runGoOra', 'context'] = "DiffPeakAnalysis"
    @params['enrichMaxPeaks'] = '2000'
    @params['enrichMaxPeaks', 'description'] = 'maximum number of top candidates per direction used for motif and GO enrichment'
    @params['enrichMaxPeaks', 'context'] = "DiffPeakAnalysis"
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
