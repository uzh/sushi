#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class HomerDiffPeaksApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'HomerDiffPeaks'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'ATAC'
    @description =<<-EOS
    Finding Peaks and Differential Peaks with or without Replicates from HOMER. <br/>
    It requires bam file mapped to UCSC reference genome (chr1 format). <br/>

    <a href='http://homer.ucsd.edu/homer/ngs/peaksReplicates.html'/>HOMER web-site</a>
EOS
    @required_columns = ['Name','BAM']
    @required_params = ['grouping', 'sampleGroup', 'refGroup']
    @params['cores'] = '4'
    @params['ram'] = '15'
    @params['scratch'] = '100'
    @params['paired'] = true
    @params['grouping'] = ''
    @params['sampleGroup'] = ''
    @params['sampleGroup', 'description'] = 'sampleGroup should be different from refGroup'
    @params['refGroup'] = ''
    @params['refGroup', 'description'] = 'refGroup should be different from sampleGroup'
    @params['refBuildHOMER'] = ['hg38', 'mm10']
    @params['refBuildHOMER', 'description'] = 'The current supported genomes from HOMER. More is available.'
    @params['repFoldChange'] = '2'
    @params['repFoldChange', 'description'] = 'Replicate fold change cutoff for peak identification (calculated by DESeq2)'
    @params['repFDR'] = '0.05'
    @params['repFDR', 'description'] = 'Replicate FDR cutoff for peak identification (calculated by DESeq2)'
    @params['balanced'] = true
    @params['balanced', 'description'] = 'Do not force the use of normalization factors to match total mapped reads.  This can be useful when analyzing differential peaks between similar data (for example H3K27ac) where we expect similar levels in all experiments. Applying this allows the data to essentially be quantile normalized during the differential calculation.'
    @params['degMethod'] = ['edgeR', 'DESeq2']
    @params['style'] = ['histone', 'factor', 'tss', 'groseq', 'dnase', 'super', 'mC']
    @params['style', 'description'] = 'Style of peaks found by findPeaks during features selection'
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'to define batches in the analysis to perform paired test, e.g. -batch 1 2 1 2'
    @params['mail'] = ""
    @modules = ["Dev/R", "Tools/HOMER", "Tools/samtools", "Tools/BEDTools"]
    @inherit_columns = ["Order Id"]
  end
  def next_dataset
    @comparison = "#{@params['sampleGroup']}--over--#{@params['refGroup']}"
    @params['comparison'] = @comparison
    @params['name'] = @comparison
    report_file = File.join(@result_dir, "#{@params['comparison']}")
    diffPeak_link = File.join(report_file, 'diffPeaks.txt')
    {'Name'=>@params['name'],
     'Report [File]'=>report_file,
     'DiffPeak [Link]'=>diffPeak_link,
    }.merge(extract_columns(colnames: @inherit_columns))
  end
  def commands
    run_RApp("EzAppHomerDiffPeaks")
  end
end

if __FILE__ == $0
  
end
