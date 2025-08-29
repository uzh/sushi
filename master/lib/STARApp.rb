#!/usr/bin/env ruby
# encoding: utf-8
Version = '20240711-164725'

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
    @params['cores'] = [8, 1, 2, 4, 8]
    @params['ram'] = ['30', '40', '60', '330']
    @params['ram', 'description'] = '30/40 GB is enough for human and mouse data, 60 GB for a hybrid genome. Only use the 200GB option for very large genomes (e.g. plants)'
    @params['scratch'] = '100'
    @params['refBuild'] = ref_selector
    @params['refBuild', 'description'] = 'Select the reference genome you wish to map your reads to. Use the most recent version by default, unless you are matching a previous run.'
    @params['paired'] = false
    @params['paired', 'description'] = 'If this is not autopopulated, check which sequencing config was used to determine. If you only have R1, set to false. If you have R1 and R2, set to true.'
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['strandMode', 'description'] = 'If this is not autopopulated, check which library kit was used to determine. If you are unsure, ask your coach.'
    @params['refFeatureFile'] = 'genes.gtf'
    @params['secondRef'] = ''
    @params['secondRef', 'description'] = 'extra DNA/RNA sequences to use for alignment; needs to point to a file on FGCZ servers; ask for upload sushi@fgcz.ethz.ch '
    @params['cmdOptions'] = '--sjdbOverhang 150 --outFilterType BySJout --outFilterMatchNmin 30 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.05 --outMultimapperOrder Random --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignIntronMax 100000 --alignMatesGapMax 100000  --outFilterMultimapNmax 50 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --outSAMstrandField intronMotif --alignEndsProtrude 10 ConcordantPair --outSAMmultNmax 4'
    @params['getJunctions'] = false
    #@params['getChimericJunctions'] = false
    @params['twopassMode'] = false
    @params['twopassMode', 'description'] = 'Per-sample 2-pass mapping or 1-pass mapping in STAR. 2-pass mapping allows to detect more splices reads mapping to novel junctions.'
    # trimming options
    # general
    #@params['trimAdapter', 'hr'] = true
    @params['trimAdapter', 'hr-header'] = "fastp parameters"
    @params['trimAdapter'] = true
    # Fastp
    ## trimming
    @params['trim_front1'] = '0'
    @params['trim_front1','description'] = 'trimming how many bases in front for read1 (and read2), default is 0.'
    @params['trim_tail1'] = '0'
    @params['trim_tail1','description'] = 'trimming how many bases in tail for read1 (and read2), default is 0.'
    @params['cut_front'] = false
    @params['cut_front','description'] = 'move a sliding window from front (5p) to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.'
    @params['cut_front_window_size'] = '4'
    @params['cut_front_mean_quality'] = '20'
    @params['cut_tail'] = false
    @params['cut_tail','description'] = 'move a sliding window from tail (3p) to front, drop the bases in the window if its mean quality < threshold, stop otherwise.'
    @params['cut_tail_window_size'] = '4'
    @params['cut_tail_mean_quality'] = '20'
    @params['cut_right'] = false
    @params['cut_right','description'] = 'move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, and then stop.'
    @params['cut_right_window_size'] = '4'
    @params['cut_right_mean_quality'] = '20'    
    @params['average_qual'] = '0'
    @params['average_qual','description'] = 'if one read average quality score <avg_qual>, then this read/pair is discarded. Default 0 means no requirement'
    @params['max_len1'] = '0'
    @params['max_len1','description'] = 'if read1 is longer than max_len1, then trim read1 at its tail to make it as long as max_len1. Default 0 means no limitation. If two reads are present, the same will apply to read2.'
    @params['max_len2'] = '0'
    @params['max_len2','description'] = 'if read1 is longer than max_len2, then trim read2 at its tail to make it as long as max_len1. Default 0 means no limitation.'
    @params['poly_x_min_len'] = '10'
    @params['poly_x_min_len','description'] = 'the minimum length to detect polyX in the read tail. 10 by default.'
    @params['length_required'] = '18'
    @params['length_required','description'] = 'reads shorter than length_required will be discarded.'
    @params['cmdOptionsFastp'] = ''
    @params['barcodePattern', 'hr-header'] = 'UMI tools'
    @params['barcodePattern'] = '' 
    @params['barcodePattern','description'] = 'optional for libraries which are including UMIs e.g. NNNNNNNN for TaKaRa SMARTer pico RNA kit v3' 
    ## additional commands
    @params['markDuplicates', 'hr-header'] = 'Additional parameters'
    @params['markDuplicates'] = false
    @params['markDuplicates', 'description'] = 'should duplicates be marked with picard. It is recommended for ChIP-seq and ATAC-seq data.'
    @params['specialOptions'] = ''

    @params['mail'] = ""
    # Python2 is required because of RSeQC package
    @modules = ["Aligner/STAR", "Tools/samtools", "Dev/jdk", "Dev/R", "Dev/Python", "Tools/Picard", "QC/fastp"]
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
        'IGV [Link,File]'=>File.join(@result_dir, "#{@dataset['Name']}-igv.html"),
        'Species'=>@dataset['Species'],
        'refBuild'=>@params['refBuild'],
        'paired'=>@params['paired'],
        'refFeatureFile'=>@params['refFeatureFile'],
        'strandMode'=>@params['strandMode'],
        'Read Count'=>@dataset['Read Count'],
        'PreprocessingLog [File]'=>File.join(@result_dir, "#{@dataset['Name']}_preprocessing.log"),
        'STARLog [File]'=>File.join(@result_dir, "#{@dataset['Name']}_STAR.log")
    }.merge(extract_columns(tags: @inherit_tags))
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
