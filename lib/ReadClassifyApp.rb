#!/usr/bin/env ruby
# encoding: utf-8
Version = '20141117-013020'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class ReadClassifyApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'ReadClassify'
    @description =<<-EOS
ReadClassify does the mapped read classification depending on the number of mismatch. ReadClassifyApp calls 'read_clasify.py' script.
Three types of bam files will be generated.<br />
Panret1Orig: the read mapped on build1 with less mismatch than build2<br />
Panret1Other: the read mapped on build1 with more mismatch than build2<br />
Panret1Common: the read mapped on build1 with the same number of mismatch to build2<br />
Panret2XX: vice versa<br />
Panret1Common and Parent2Common should be the same but Parent1Orig and Panret2Other are usually different.<br />
http://seselab.org/homeoroq/
    EOS

    @analysis_category = 'Demo'
    @required_columns = ['Name','BAM1','BAM2','build1','build2','Species']
    @required_params = []
    # optional params
    @params['cores'] = '1'
    @params['ram'] = '16'
    @params['scratch'] = '100'
  end
  def next_dataset
    {'Name'=>@dataset['Name'], 
     'Parent1Orig [File]'=>File.join(@result_dir, "parent1_genome_orig.bam"), 
     'Parent1Other [File]'=>File.join(@result_dir, "parent1_genome_other.bam"), 
     'Parent1Common [File]'=>File.join(@result_dir, "parent1_genome_common.bam"), 
     'Parent2Orig [File]'=>File.join(@result_dir, "parent2_genome_orig.bam"), 
     'Parent2Other [File]'=>File.join(@result_dir, "parent2_genome_other.bam"), 
     'Parent2Common [File]'=>File.join(@result_dir, "parent2_genome_common.bam"), 
     'Species'=>@dataset['Species'],
    }.merge factor_dataset
  end
  def commands
    bam1 = File.join(@gstore_dir, @dataset['BAM1'])
    bam2 = File.join(@gstore_dir, @dataset['BAM2'])
    sam1 = File.join("tmp", File.basename(bam1).gsub(/.bam/, '_parent1.sam'))
    sam2 = File.join("tmp", File.basename(bam2).gsub(/.bam/, '_parent2.sam'))
    out1_prefix = "parent1_genome"
    out2_prefix = "parent2_genome"
    command = "mkdir tmp\n"
    command << "samtools view -h #{bam1} > #{sam1}\n"
    command << "samtools view -h #{bam2} > #{sam2}\n"
    command << "read_classify.py #{sam1} #{sam2} #{out1_prefix} #{out2_prefix}\n"
    command
  end
end


