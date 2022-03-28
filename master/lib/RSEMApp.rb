#!/usr/bin/env ruby
# encoding: utf-8
Version = '20180530-104624'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class RSEMApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'RSEM'
    @analysis_category = 'Count'
    @description =<<-EOS
    Use bowtie alignments to transcript database and a posterior model to estimate isoform/gene abundances<br/>
<a href='http://deweylab.biostat.wisc.edu/rsem/rsem-calculate-expression.html'>manual/</a><br/>
Noteworthy is the option --bowtie-e which can be used to limit the sum of mismatching qualities for the alignments
EOS
    @required_columns = ['Name','Read1','Species']
    @required_params = ['paired', 'strandMode']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '100'
    @params['refBuild'] = ref_selector
    @params['paired'] = false
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['refFeatureFile'] = 'genes.gtf'
    @params['bowtie-e'] = '200'
    @params['bowtie-e', 'description'] = 'maximum sum of base qualities at mismatching positions'
    @params['cmdOptions'] = ' --calc-ci --sort-bam-by-read-name'
    @params['keepBam'] = false
    @params['keepBam', 'description'] = 'converts the transcript alignments into genome coordinates and reports them as a BAM file'
    
    
    @params['transcriptFasta'] = ''
    @params['transcriptFasta', 'description'] = 'give full path of transcript fasta file; in that case the build is ignored; if it comes from trinity assembly the gene-isoform associations will be extracted and used'
    @params['transcriptTypes'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA', 'long_noncoding', 'short_noncoding', 'pseudogene']
    @params['transcriptTypes', 'multi_selection'] = true
    @params['transcriptTypes', 'selected'] = 'protein_coding'
    
    # trimming options
    # general
    @params['trimAdapter', 'hr'] = true
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
    ## additional commands
    @params['markDuplicates'] = true
    @params['markDuplicates', 'description'] = 'should duplicates be marked with picard. It is recommended for ChIP-seq and ATAC-seq data.'
    
    @params['mail'] = ""
	  # Bowtie >=1.2.0 may return interleaving mates which trips RSEM (as of v1.3.0) as it expects
	  # each read to be followed by a mate.
	  # Also, as of v1.3.0, it only supports samtools v1.3.1
    @modules = ["Tools/samtools", "Aligner/Bowtie", "Aligner/Bowtie2", "Aligner/STAR", "Aligner/RSEM", "QC/fastp", "Dev/R"]
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
    dataset = {
      'Name'=>@dataset['Name'],
      'Count [File]'=>File.join(@result_dir, "#{@dataset['Name']}.txt")
    }
    keep_bam = {
      'BAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam"),
      'BAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam.bai")
    }
    non_keep_bam = {
      'PreprocessingLog [File]'=>File.join(@result_dir, "#{@dataset['Name']}_preprocessing.log")
    }
    non_transcript_fasta = {
      'transcriptTypes'=>@params['transcriptTypes']
    }
    transcript_fasta = {
      'transcriptTypes'=>'NA'
    }
    common = {
      'Species'=>@dataset['Species'],
      'refBuild'=>@params['refBuild'],
      'featureLevel'=>'isoform',
      'refFeatureFile'=>@params['refFeatureFile'],
      'strandMode'=>@params['strandMode'],
      'paired'=>@params['paired'],
      'Read Count'=>@dataset['Read Count'],
    }
    merge_transcript_fasta =->() do
      if @params['transcriptFasta'].to_s.empty?
        dataset.merge!(non_transcript_fasta)
      else
        dataset.merge!(transcript_fasta)
      end
    end
    if @params['keepBam']
      dataset.merge!(keep_bam)
      dataset.merge!(common)
      merge_transcript_fasta.()
    else
      dataset.merge!(common)
      merge_transcript_fasta.()
      dataset.merge!(non_keep_bam)
    end
    dataset.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppRSEM")
  end
end

if __FILE__ == $0
  run RSEMApp

end
