#!/usr/bin/env ruby
# encoding: utf-8
Version = '20190709-112022'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class EAGLERCApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'EAGLERC'
    @description =<<-EOS
EAGLE: Explicit Alternative Genome Likelihood Evaluator <br />
<br />
* https://github.com/tony-kuo/eagle
    EOS

    @analysis_category = 'Polyploid'
    @required_columns = ['Name','BAM1','BAM2','refBuild1','refBuild2','Species']
    @required_params = []
    # optional params
    @params['cores'] = '1'
    @params['ram'] = '16'
    @params['scratch'] = '100'
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
    @modules = ["Tools/EAGLE/1.1.1"]
  end
  def preprocess
    # first, from species name
    @parent1_genome = if samp1 = @dataset_hash.first and ref_path1 = samp1['refBuild1'] and dirs1 = ref_path1.split('/') and spc1 = dirs1.first and sub1 = spc1.split('_')
                        spc1[0] + sub1.last[0,3]
                      else
                        '1'
                      end
    @parent2_genome = if samp2 = @dataset_hash.first and ref_path2 = samp2['refBuild2'] and dirs2 = ref_path2.split('/') and spc2 = dirs2.first and sub2 = spc2.split('_')
                        spc2[0] + sub2.last[0,3]
                      else
                        '2'
                      end
    if @parent1_genome == @parent2_genome
      # second, from build name
      build1_parts = if build1 = dirs1[2]
                       build1.split('_')
                     else
                       ["1"]
                     end
      build2_parts = if build2 = dirs2[2]
                       build2.split('_')
                     else
                       ["2"]
                     end
      @parent1_genome = build1_parts.join("_")
      @parent2_genome = build2_parts.join("_")
    end
    if @parent1_genome == @parent2_genome
      # otherwise, 1 and 2
      @parent1_genome = "1"
      @parent2_genome = "2"
    end
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'Parent1RefBAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}_#{@parent1_genome}_ref.bam"),
     'Parent1AltBAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}_#{@parent1_genome}_alt.bam"),
     'Parent1UnkBAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}_#{@parent1_genome}_unk.bam"),
     'Parent1MulBAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}_#{@parent1_genome}_mul.bam"),
     'refBuild1'=>@dataset['refBuild1'],
     'Parent2RefBAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}_#{@parent2_genome}_ref.bam"),
     'Parent2AltBAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}_#{@parent2_genome}_alt.bam"),
     'Parent2UnkBAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}_#{@parent2_genome}_unk.bam"),
     'Parent2MulBAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}_#{@parent2_genome}_mul.bam"),
     'stdout log [File]'=>File.join(@result_dir, "#{@dataset['Name']}.sort.stdout.log"),
     'errout log [File]'=>File.join(@result_dir, "#{@dataset['Name']}.sort.errout.log"),
     'refBuild2'=>@dataset['refBuild2'],
     'Species'=>@dataset['Species'],
     'dummy [File]'=>File.join(@result_dir, "#{@dataset['Name']}_dummy.txt")
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    bam1 = File.join(@gstore_dir, @dataset['BAM1'])
    bam2 = File.join(@gstore_dir, @dataset['BAM2'])

    command = "eagle --version\n"
    ref1 = File.join(GENOME_REF_DIR, @dataset['refBuild1'].split('/')[0,3].join("/")+"/Sequence/WholeGenomeFasta/genome.fa")
    ref2 = File.join(GENOME_REF_DIR, @dataset['refBuild2'].split('/')[0,3].join("/")+"/Sequence/WholeGenomeFasta/genome.fa")
    command << "eagle-rc --ngi --ref1=#{ref1} --ref2=#{ref2} --bam1=#{bam1} --bam2=#{bam2} -o #{@dataset['Name']} > #{@dataset['Name']}.sort.stdout.log 2> #{@dataset['Name']}.sort.errout.log\n"
    command << "mv #{@dataset['Name']}1.ref.bam #{@dataset['Name']}_#{@parent1_genome}_ref.bam\n"
    command << "mv #{@dataset['Name']}1.alt.bam #{@dataset['Name']}_#{@parent1_genome}_alt.bam\n"
    command << "mv #{@dataset['Name']}1.unk.bam #{@dataset['Name']}_#{@parent1_genome}_unk.bam\n"
    command << "mv #{@dataset['Name']}1.mul.bam #{@dataset['Name']}_#{@parent1_genome}_mul.bam\n"

    command << "mv #{@dataset['Name']}2.ref.bam #{@dataset['Name']}_#{@parent2_genome}_ref.bam\n"
    command << "mv #{@dataset['Name']}2.alt.bam #{@dataset['Name']}_#{@parent2_genome}_alt.bam\n"
    command << "mv #{@dataset['Name']}2.unk.bam #{@dataset['Name']}_#{@parent2_genome}_unk.bam\n"
    command << "mv #{@dataset['Name']}2.mul.bam #{@dataset['Name']}_#{@parent2_genome}_mul.bam\n"

    command << "echo '#{GlobalVariables::SUSHI}' > #{@dataset['Name']}_dummy.txt\n"
    command
  end
end


