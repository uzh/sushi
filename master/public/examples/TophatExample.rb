#!/usr/bin/env ruby
# encoding: utf-8
Version = '20160426-053404'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class TophatExample < SushiFabric::SushiApp
  def initialize
    super
    @name = 'TophatExample'
    @analysis_category = 'Map'
    @description =<<-EOS
TopHat is a fast splice junction mapper for RNA-Seq reads. 
It aligns RNA-Seq using bowtie2, and then analyzes the mapping results to identify splice junctions between exons.<br />
<a href='https://ccb.jhu.edu/software/tophat/manual.shtml'>https://ccb.jhu.edu/software/tophat/manual.shtml/</a>
    EOS
    @required_columns = ['Name','Read1','Species']
    @required_params = ['build','paired']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '15'
    @params['scratch'] = '100'
    @params['library_type'] = ['fr-unstranded', 'fr-firststrand', 'fr-secondstrand']
    @params['library_type', 'description'] = 'strand specificity of the library'
    @params['paired'] = false
    @params['paired', 'description'] = 'either the reads are paired-end or single-end'
    @params['build'] = ['GRCh38', 'hg19', 'GRCm38', 'mm10']
    @params['build', 'description'] = 'build id of the genome assembly'
    @params['tophat_options'] = ""
    @params['tophat_options', 'description'] = "other options that will be directly pasted to tophats command line"
  end
  def set_default_parameters
    @params['paired'] = dataset_has_column?('Read2')
  end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  end
  def next_dataset
    {'Name'=>@dataset['Name'], 
     'BAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam"), 
     'BAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam.bai"),
     'build'=>@params['build'],
     'Species'=>@dataset['Species']
    }
  end
  def bowtie2_index
    @bowtie2_index= File.join('/srv/GT/reference', @params['build'], 'BOWTIE2Index/genome')
  end
  def num_threads
    if @params['cores'].to_i > 1
      "--num-threads #{@params['cores']}"
    else
      ""
    end 
  end
  def commands
    cmd = "which tophat; tophat -v\n"
    cmd << "tophat -o . #{num_threads} --library-type #{@params['library_type']} #{@params['tophat_options']} #{bowtie2_index} #{@gstore_dir}/#{@dataset['Read1']}"
    if @params['paired']
      cmd << ",#{@gstore_dir}/#{@dataset['Read2']}"
    end
    cmd << "\n"
    cmd << "mv accepted_hits.bam #{@dataset['Name']}.bam\n"
    cmd << "samtools index #{@dataset['Name']}.bam\n"
  end
end
