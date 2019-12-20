#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class ExceRptApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Excerpt'
    @params['process_mode'] = 'SAMPLE'
    @analysis_category = 'Count'
    @description =<<-EOS
Annotation and Profiling of smallRNA-seq with exceRpt's pipeline<br/>
ExceRpt performs a complete analysis of small-RNA-seq. This software carries out processing, filtering, and alignment for each smallRNA-seq sample from a run.<br />
<a href='https://github.com/rkitchen/exceRpt'>https://github.com/rkitchen/exceRpt/</a>
EOS
    @required_columns = ['Name','Read1', 'Adapter1', 'Species']
    @required_params = ['refBuild', 'cores', 'ram', 'scratch']
    
    # parameters
    ## general
    @params['cores'] = '8'
    @params['ram'] = '16'
    #@params['JAVA_RAM'] = '10G'
    @params['scratch'] = '100'
    @params['refBuild'] = ['Homo_sapiens/UCSC/hg38', 'Mus_musculus/UCSC/mm10']  ## this param substitutes exceRpt_smallRNA's vaibale: "MAIN_ORGANISM_GENOME_ID"
    @params['mail'] = ""
    
    ## tunable parameters exceRpt_smallRNA
    #@params['ADAPTER_SEQ'] = 'guessKnown' # taken from dataset
    @params['SAMPLE_NAME'] = 'NULL'
    @params['SAMPLE_NAME','description'] = 'If "NULL", the ExceRptApp will automatically set the name from the fastq(.gz) files.'
    @params['ENDOGENOUS_LIB_PRIORITY'] = 'miRNA,tRNA,piRNA,gencode,circRNA'
    @params['ENDOGENOUS_LIB_PRIORITY','description'] = '<comma,separated,list,no,spaces> to choose the priority of each library during read assignment and quantification.'
    ## preprocessing
    @params['TRIM_N_BASES_5p'] = '0'
    @params['TRIM_N_BASES_3p'] = '0'
    @params['MIN_READ_LENGTH'] = '18'
    @params['QFILTER_MIN_QUAL'] = '20'
    @params['QFILTER_MIN_READ_FRAC'] = '80'
    ## barcodes
    @params['RANDOM_BARCODE_LENGTH'] = '0'
    @params['RANDOM_BARCODE_LOCATION'] = ['-5p -3p','-5p','-3p']
    @params['KEEP_RANDOM_BARCODE_STATS'] = ['false','true']
    ## calibration (to be further tested before making it available)
    #@params['CALIBRATOR_LIBRARY'] = 'NULL'
    #@params['CALIBRATOR_LIBRARY','description'] = 'Path to a bowtie2 index of calibrator oligos spiked-in for QC or normalisation.'
    #@params['CALIBRATOR_TRIM5p'] = '0'
    #@params['CALIBRATOR_TRIM3p'] = '0'
    
    @params['DOWNSAMPLE_RNA_READS'] = 'NULL'
    @params['DOWNSAMPLE_RNA_READS','description'] = 'Choose whether to downsample to this number of reads after assigning reads to the various transcriptome libraries (may be useful for normalising very different yields)'
    #@params['MAP_EXOGENOUS'] = ['off','miRNA'] # we do not use this functionality
    
    ## STAR parameters
    @params['STAR_alignEndsType'] = ['Local','EndToEnd']
    @params['STAR_alignEndsType','description'] = 'Defines the alignment mode; local alignment is recommended to allow for isomiRs.'
    @params['STAR_outFilterMatchNmin'] = '18'
    @params['STAR_outFilterMatchNmin','description'] = 'Minimum number of bases to include in the alignment (should match the minimum read length defined above).'
    @params['STAR_outFilterMatchNminOverLread'] = '0.9'
    @params['STAR_outFilterMatchNminOverLread','description'] = 'Minimum fraction of the read that *must* remain following soft-clipping in a local alignment.'
    @params['STAR_outFilterMismatchNmax'] = '1'
    @params['STAR_outFilterMismatchNmax','description'] = 'Maximum allowed mismatched bases in the aligned portion of the read.'
    #@params['MAX_MISMATCHES_EXOGENOUS'] = '0'
    @params['REMOVE_LARGE_INTERMEDIATE_FILES'] = ['true','false']
    @params['REMOVE_LARGE_INTERMEDIATE_FILES','description'] = 'When exceRpt finishes, choose whether to remove the large alignment files that can take a lot of disk space.'

    ## modules
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
    @modules = ["Dev/R"]
  end
  def next_dataset
    dataset = {
      'Name'=>@dataset['Name'],
      'excerpt [File]'=>File.join(@result_dir, "#{@dataset['Name']}"),
      'Species'=>@dataset['Species'],
      'refBuild'=>@params['refBuild']
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppExceRpt")
  end
end

if __FILE__ == $0
  run ExceRptApp
end
