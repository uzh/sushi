#!/usr/bin/env ruby
# encoding: utf-8
Version = '20150524-210956'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class ReadClassifyApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'ReadClassify'
    @description =<<-EOS
ReadClassifyApp does the mapped read classification depending on the number of mismatch. ReadClassifyApp calls 'read_clasify.py' script.
Three types of bam files will be generated.<br />
Panret1Orig: the read mapped on refBuild1 with less mismatch than refBuild2<br />
Panret1Other: the read mapped on refBuild1 with more mismatch than refBuild2<br />
Panret1Common: the read mapped on refBuild1 with the same number of mismatch to refBuild2<br />
Panret2XX: vice versa<br />
Panret1Common and Parent2Common should be the same but Parent1Orig and Panret2Other are usually different.<br />
http://seselab.org/homeoroq/
    EOS

    @analysis_category = 'Demo'
    @required_columns = ['Name','BAM1','BAM2','refBuild1','refBuild2','Species']
    @required_params = []
    # optional params
    @params['cores'] = '1'
    @params['ram'] = '16'
    @params['scratch'] = '100'
  end
  def next_dataset
    {'Name'=>@dataset['Name'], 
     'Parent1OrigBAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}_parent1_genome_orig.bam"), 
     'Parent1OrigBAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}_parent1_genome_orig.bam.bai"), 
     'Parent1OtherBAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}_parent1_genome_other.bam"), 
     'Parent1OtherBAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}_parent1_genome_other.bam.bai"), 
     'Parent1CommonBAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}_parent1_genome_common.bam"), 
     'Parent1CommonBAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}_parent1_genome_common.bam.bai"), 
     'build1'=>@dataset['build1'],
     'Parent2OrigBAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}_parent2_genome_orig.bam"), 
     'Parent2OrigBAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}_parent2_genome_orig.bam.bai"), 
     'Parent2OtherBAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}_parent2_genome_other.bam"), 
     'Parent2OtherBAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}_parent2_genome_other.bam.bai"), 
     'Parent2CommonBAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}_parent2_genome_common.bam"), 
     'Parent2CommonBAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}_parent2_genome_common.bam.bai"), 
     'build2'=>@dataset['build2'],
     'Species'=>@dataset['Species'],
    }.merge(extract_column("Factor")).merge(extract_column("B-Fabric"))
  end
  def commands
    bam1 = File.join(@gstore_dir, @dataset['BAM1'])
    bam2 = File.join(@gstore_dir, @dataset['BAM2'])
    sam1 = File.join("tmp", File.basename(bam1).gsub(/.bam/, '_parent1.sam'))
    sam2 = File.join("tmp", File.basename(bam2).gsub(/.bam/, '_parent2.sam'))
    out1_prefix = "#{@dataset['Name']}_parent1_genome"
    out2_prefix = "#{@dataset['Name']}_parent2_genome"
    command = "mkdir tmp\n"
    command << "grep 'version' /usr/local/ngseq/stow/read_classify-2.1.0/bin/read_classify.py\n"
    command << "/usr/local/ngseq/bin/python --version\n"
    command << "samtools view -h #{bam1} > #{sam1}\n"
    command << "samtools view -h #{bam2} > #{sam2}\n"
    command << "/usr/local/ngseq/bin/python /usr/local/ngseq/stow/read_classify-2.1.0/bin/read_classify.py #{sam1} #{sam2} #{out1_prefix} #{out2_prefix}\n"
    ["orig", "other", "common"].each do |type|
      command << "samtools index #{out1_prefix}_#{type}.bam\n" 
      command << "samtools index #{out2_prefix}_#{type}.bam\n" 
    end
    command
  end
end


