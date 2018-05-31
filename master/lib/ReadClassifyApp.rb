#!/usr/bin/env ruby
# encoding: utf-8
Version = '20180531-092244'

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
The reads included in Panret1Common and Parent2Common should be same but the reference information is different in the BAM files<br />
Parent1Orig and Parent2Other, Parent1Other and Parent2Orig, are as well<br />
http://seselab.org/homeoroq/
    EOS

    @analysis_category = 'Polyploid'
    @required_columns = ['Name','BAM1','BAM2','refBuild1','refBuild2','Species']
    @required_params = []
    # optional params
    @params['cores'] = '1'
    @params['ram'] = '16'
    @params['scratch'] = '100'
    @params['paired'] = true
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
    @modules = ["Dev/Python/2.7.13", "Tools/samtools", "Tools/sambamba"]
  end
  def preprocess
    @parent1_genome = if samp = @dataset_hash.first and ref_path = samp['refBuild1'] and dirs = ref_path.split('/') and spc = dirs.first and sub = spc.split('_')
                spc[0] + sub.last[0,3]
              else
                'parent1_genome'
              end
    @parent2_genome = if samp = @dataset_hash.first and ref_path = samp['refBuild2'] and dirs = ref_path.split('/') and spc = dirs.first and sub = spc.split('_')
                     spc[0] + sub.last[0,3]
                   else
                     'parent2_genome'
                   end
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'Parent1OrigBAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}_#{@parent1_genome}_orig.bam"),
     'Parent1OrigBAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}_#{@parent1_genome}_orig.bam.bai"),
     'Parent1OtherBAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}_#{@parent1_genome}_other.bam"),
     'Parent1OtherBAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}_#{@parent1_genome}_other.bam.bai"),
     'Parent1CommonBAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}_#{@parent1_genome}_common.bam"),
     'Parent1CommonBAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}_#{@parent1_genome}_common.bam.bai"),
     'refBuild1'=>@dataset['refBuild1'],
     'Parent2OrigBAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}_#{@parent2_genome}_orig.bam"),
     'Parent2OrigBAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}_#{@parent2_genome}_orig.bam.bai"),
     'Parent2OtherBAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}_#{@parent2_genome}_other.bam"),
     'Parent2OtherBAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}_#{@parent2_genome}_other.bam.bai"),
     'Parent2CommonBAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}_#{@parent2_genome}_common.bam"),
     'Parent2CommonBAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}_#{@parent2_genome}_common.bam.bai"),
     'refBuild2'=>@dataset['refBuild2'],
     'Species'=>@dataset['Species'],
     'dummy [File]'=>File.join(@result_dir, "#{@dataset['Name']}_dummy.txt")
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    bam1 = File.join(@gstore_dir, @dataset['BAM1'])
    bam2 = File.join(@gstore_dir, @dataset['BAM2'])
    sam1 = File.join("tmp", File.basename(bam1).gsub(/.bam/, '_parent1.sam'))
    sam2 = File.join("tmp", File.basename(bam2).gsub(/.bam/, '_parent2.sam'))
    out1_prefix = "#{@dataset['Name']}_#{@parent1_genome}"
    out2_prefix = "#{@dataset['Name']}_#{@parent2_genome}"
    if out1_prefix == out2_prefix
      command << "echo 'ERROR: the top 3 chars of subspecies name of two parental reference folder should be different'"
      command << "exit"
    end
    command = "mkdir tmp\n"
    command << "grep 'version' /usr/local/ngseq/stow/read_classify-2.1.0/bin/read_classify.py\n"
    command << "python --version\n"
    #command << "export PYTHONPATH=$PYTHONPATH:/usr/local/ngseq/lib/python2.7:/usr/local/ngseq/lib/python2.7/dist-packages\n"
    command << "env|grep PYTHON\n"
    command << "samtools view -h #{bam1} > #{sam1}\n"
    command << "samtools view -h #{bam2} > #{sam2}\n"

    if @params['paired']
      command << "python /usr/local/ngseq/stow/read_classify-2.1.0/bin/read_classify.py #{sam1} #{sam2} #{out1_prefix} #{out2_prefix}\n"
    else
      command << "python /usr/local/ngseq/stow/read_classify-2.1.0/bin/read_classify_single_end.py #{sam1} #{sam2} #{out1_prefix} #{out2_prefix}\n"
    end
    ["orig", "other", "common"].each do |type|
      command << "samtools index #{out1_prefix}_#{type}.bam\n"
      command << "samtools index #{out2_prefix}_#{type}.bam\n"
    end
    command << "echo '#{GlobalVariables::SUSHI}' > #{@dataset['Name']}_dummy.txt\n"
    command
  end
end


