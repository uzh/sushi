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
<a href='https://github.com/rkitchen/exceRpt'>https://rkitchen.github.io/exceRpt/</a>
EOS
    @required_columns = ['Name','Read1', 'Adapter1', 'Species']
    @required_params = ['refBuild', 'cores', 'ram', 'scratch']
    
    # parameters
    ## general
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '100'
    @params['refBuild'] = ['Homo_sapiens/UCSC/hg38', 'Mus_musculus/UCSC/mm10']  ## this param substitutes exceRpt_smallRNA's vaibale: "MAIN_ORGANISM_GENOME_ID"

    ## tunable parameters exceRpt_smallRNA
    @params['ENDOGENOUS_LIB_PRIORITY'] = 'miRNA,tRNA,piRNA,gencode,circRNA'
    @params['ENDOGENOUS_LIB_PRIORITY','description'] = 'Choose the priority of each library during read assignment and quantification, e.g. "comma,separated,list,with,no,spaces".'
    
    ## preprocessing
    @params['TRIM_N_BASES_5p'] = '0'
    @params['TRIM_N_BASES_3p'] = '0'
    @params['MIN_READ_LENGTH'] = '18'
    @params['QFILTER_MIN_QUAL'] = '20'
    @params['QFILTER_MIN_READ_FRAC'] = '80'
    
    ## barcodes
    @params['RANDOM_BARCODE_LENGTH'] = '0'
    @params['RANDOM_BARCODE_LENGTH','description'] = 'Identify and remove random barcodes of this number of nucleotides. For a Bioo prep with a 4N random barcode on both the 3p and 5p adapter, this value should be 4.'
    @params['RANDOM_BARCODE_LOCATION'] = ['-5p -3p','-5p','-3p']
    @params['RANDOM_BARCODE_LOCATION','description'] = 'Specify where to look for the random barcode(s).'
    @params['KEEP_RANDOM_BARCODE_STATS'] = false
    @params['KEEP_RANDOM_BARCODE_STATS', 'description'] = 'Specify whether or not to calculate overrepresentation statistics using the random barcodes (this may be slow and memory intensive!).'
    
    ## STAR parameters
    @params['STAR_alignEndsType'] = ['Local','EndToEnd']
    @params['STAR_alignEndsType','description'] = 'Define the alignment mode; local alignment is recommended to allow for isomiRs.'
    @params['STAR_outFilterMatchNmin'] = '18'
    @params['STAR_outFilterMatchNmin','description'] = 'Minimum number of bases to include in the alignment (should match the minimum read length defined above).'
    @params['STAR_outFilterMatchNminOverLread'] = '0.9'
    @params['STAR_outFilterMatchNminOverLread','description'] = 'Minimum fraction of the read that *must* remain following soft-clipping in a local alignment.'
    @params['STAR_outFilterMismatchNmax'] = '1'
    @params['STAR_outFilterMismatchNmax','description'] = 'Maximum allowed mismatched bases in the aligned portion of the read.'
    @params['REMOVE_LARGE_INTERMEDIATE_FILES'] = ['true','false']
    @params['REMOVE_LARGE_INTERMEDIATE_FILES','description'] = 'When exceRpt finishes, choose whether to remove the large alignment files that can take a lot of disk space.'
    
    @params['cmdOptions'] = 'JAVA_RAM=10G'

    @params['mail'] = ""

    ## modules
    @params['ExcerptVersion'] = ["Tools/exceRpt/5.0.0", "Tools/exceRpt/4.6.3"]
    @params['starVersion'] = ["Aligner/STAR/2.6.1e"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
    @modules = ["Dev/R", "Aligner/Bowtie2", "Tools/samtools", "Dev/jdk", "QC/FastQC", "Tools/fastx"]
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
    command = "module load #{@params["ExcerptVersion"]}\n"
    command << "module load #{@params["starVersion"]}\n"
    command << run_RApp("EzAppExceRpt")
  end
end

if __FILE__ == $0
  run ExceRptApp
end
