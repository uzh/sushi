#!/usr/bin/env ruby
# encoding: utf-8
Version = '20200225-143601'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class STARApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'STAR'
    @analysis_category = 'Map'
    @description =<<-EOS
    Ultafast spliced alignment<br/>
<a href='https://code.google.com/p/rna-star/'>manual</a><br/>
Noteworthy options:<ul>
<li> --outFilterMatchNmin 30 --outFilterMismatchNmax 5 --outFilterMismatchNoverLmax 0.05 --outFilterMultimapNmax 50 --alignEndsProtrude 3 ConcordantPair</li>
<li> large numbers of contigs ask form more RAM. In that case the index must be built with a smaller --genomeChrBinNbits 18;
e.g. the Wheat assembly with 1Mio contigs and --genomeChrBinNbits 10 still takes 120GB of RAM during alignmnet</li>
</ul>
EOS
    @required_columns = ['Name','Read1','Species']
    @required_params = ['refBuild','paired', 'strandMode']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '100'
    @params['refBuild'] = ref_selector
    @params['paired'] = false
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['refFeatureFile'] = 'genes.gtf'
    @params['cmdOptions'] = '--outFilterType BySJout --outFilterMatchNmin 30 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.05 --outMultimapperOrder Random --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignIntronMax 100000 --alignMatesGapMax 100000  --outFilterMultimapNmax 50 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --outSAMstrandField intronMotif --alignEndsProtrude 3 ConcordantPair'
    @params['getJunctions'] = false
    #@params['getChimericJunctions'] = false
    @params['twopassMode'] = false
    @params['twopassMode', 'description'] = 'Per-sample 2-pass mapping or 1-pass mapping in STAR. 2-pass mapping allows to detect more splices reads mapping to novel junctions.'
    # trimming options
    # general
    @params['trimAdapter', 'hr'] = true
    @params['trimAdapter'] = true
    # Fastp
    ## trimming
    @params['trim_front'] = '0'
    @params['trim_front','description'] = 'trimming how many bases in front for read1 (and read2), default is 0.'
    @params['trim_tail'] = '0'
    @params['trim_tail','description'] = 'trimming how many bases in tail for read1 (and read2), default is 0.'
    @params['cut_front'] = '0'
    @params['cut_front','description'] = 'move a sliding window from front (5p) to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.'
    @params['cut_tail'] = '0'
    @params['cut_tail','description'] = 'move a sliding window from tail (3p) to front, drop the bases in the window if its mean quality < threshold, stop otherwise.'
    @params['cut_right'] = '0'
    @params['cut_right','description'] = 'move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, and then stop.'
    @params['average_qual'] = '0'
    @params['average_qual','description'] = 'if one read average quality score <avg_qual>, then this read/pair is discarded. Default 0 means no requirement'
    @params['max_len1'] = '0'
    @params['max_len1','description'] = 'if read1 is longer than max_len1, then trim read1 at its tail to make it as long as max_len1. Default 0 means no limitation. If two reads are present, the same will apply to read2.'
    @params['poly_x_min_len'] = '10'
    @params['poly_x_min_len','description'] = 'the minimum length to detect polyX in the read tail. 10 by default.'
    @params['length_required'] = '18'
    @params['length_required','description'] = 'reads shorter than length_required will be discarded.'
    ## additional commands
    @params['cmdOptionsFastp'] = ''
    @params['markDuplicates'] = true
    @params['markDuplicates', 'description'] = 'should duplicates be marked with sambamba. It is recommended for ChIP-seq and ATAC-seq data.'
    @params['mail'] = ""
    # Python2 is required because of RSeQC package
    @modules = ["Aligner/STAR", "Tools/samtools", "QC/Flexbar", "Dev/jdk", "Tools/Picard", "QC/Trimmomatic", "Dev/Python", "Dev/R", "Tools/sambamba"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  end
  def set_default_parameters
    if dataset_has_column?('refBuild')
      @params['refBuild'] = @dataset[0]['refBuild']
    end
    @params['paired'] = dataset_has_column?('Read2')
    if dataset_has_column?('strandMode')
      @params['strandMode'] = @dataset[0]['strandMode']
    end
  end
  def next_dataset
     dataset = {
        'Name'=>@dataset['Name'],
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
        'PreprocessingLog [File]'=>File.join(@result_dir, "#{@dataset['Name']}_preprocessing.log"),
        'STARLog [File]'=>File.join(@result_dir, "#{@dataset['Name']}_STAR.log")
    }.merge(extract_columns(@inherit_tags))
     if @params['getJunctions']
       dataset['Junctions [File]'] = File.join(@result_dir, "#{@dataset['Name']}_SJ.out.tab")
       dataset['Chimerics [File]'] = File.join(@result_dir, "#{@dataset['Name']}.chimeric")
     end
     dataset
  end
  def commands
    run_RApp("EzAppSTAR")
  end
end

if __FILE__ == $0
  run STARApp

end
