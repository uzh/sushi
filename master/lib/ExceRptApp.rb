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
    ## parameters
    ### general
    @params['cores'] = '8'
    @params['ram'] = '40'
    @params['ram','description'] = 'in Gigabytes'
    #@params['JAVA_RAM'] = '10G'
    @params['scratch'] = '100'
    @params['refBuild'] = ['Mus_musculus/UCSC/mm10', 'Homo_sapiens/UCSC/hg38']
    @params['mail'] = ""
    
    ### tunable parameters exceRpt_smallRNA
    #@params['ADAPTER_SEQ'] = 'guessKnown'
    #@params['ADAPTER_SEQ','description'] = 'Write your adapter sequence or "none" in case adapters were trimmed. By default, "guessKnown" will attempt to guess the 3p adapter using known sequences.'
    @params['SAMPLE_NAME'] = 'NULL'
    @params['SAMPLE_NAME','description'] = 'If "NULL", the ExceRptApp will automatically set the name from the fastq(.gz) files.'
    #@params['MAIN_ORGANISM_GENOME_ID'] = ['hg38','mm10']
    @params['ENDOGENOUS_LIB_PRIORITY'] = 'miRNA,tRNA,piRNA,gencode,circRNA'
    @params['ENDOGENOUS_LIB_PRIORITY','description'] = '<comma,separated,list,no,spaces>'
    
    @params['CALIBRATOR_LIBRARY'] = 'NULL'
    @params['CALIBRATOR_LIBRARY','description'] = 'Path to a bowtie2 index of calibrator oligos used for QC or normalisation.'
    @params['CALIBRATOR_TRIM5p'] = '0'
    @params['CALIBRATOR_TRIM3p'] = '0'
    
    @params['TRIM_N_BASES_5p'] = '0'
    @params['TRIM_N_BASES_3p'] = '0'
    @params['RANDOM_BARCODE_LENGTH'] = '0'
    @params['RANDOM_BARCODE_LOCATION'] = ["'-5p -3p'","'-5p'","'-3p'"]
    @params['KEEP_RANDOM_BARCODE_STATS'] = ['false','true']
    @params['DOWNSAMPLE_RNA_READS'] = 'NULL'
    #@params['MAP_EXOGENOUS'] = ['off','miRNA']
    @params['MIN_READ_LENGTH'] = '18'
    @params['QFILTER_MIN_QUAL'] = '20'
    @params['QFILTER_MIN_READ_FRAC'] = '80'
    @params['STAR_alignEndsType'] = ['Local','EndToEnd']
    @params['STAR_alignEndsType','description'] = 'Defines the alignment mode; local alignment is recommended to allow for isomiRs.'
    @params['STAR_outFilterMatchNmin'] = '18'
    @params['STAR_outFilterMatchNmin','description'] = 'Minimum number of bases to include in the alignment (should match the minimum read length defined above).'
    @params['STAR_outFilterMatchNminOverLread'] = '0.9'
    @params['STAR_outFilterMatchNminOverLread','description'] = 'Minimum fraction of the read that *must* remain following soft-clipping in a local alignment.'
    @params['STAR_outFilterMismatchNmax'] = '1'
    @params['STAR_outFilterMismatchNmax','description'] = 'Maximum allowed mismatched bases in the aligned portion of the read.'
    #@params['MAX_MISMATCHES_EXOGENOUS'] = '0'


    REMOVE_LARGE_INTERMEDIATE_FILES

    ## modules
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
    @modules = ["Dev/R"]
  end
  def preprocess
    if @params['paired']
      @required_columns<<  'Read2'
    end
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
