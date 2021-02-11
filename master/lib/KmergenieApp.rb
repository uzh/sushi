#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-094757'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class KmergenieApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Kmergenie'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'QC'
    @description =<<-EOS
Kmergenie calculates kmer distribution, and it estimates the best kmer size and genome size<br/>
<a href='http://kmergenie.bx.psu.edu/'>http://kmergenie.bx.psu.edu/</a>
    EOS
    @required_columns = ['Name','Read1']
    @required_params = ['model','paired']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '15'
    @params['scratch'] = '100'
    @params['model'] = ['diploid', 'haploid']
    @params['paired'] = false
    @params['sampling'] = 'all'
    @params['sampling', 'description'] = 'all: use all reads, number: use the number of random sampled reads'
    @params['cmdOptions'] = ''
	  @modules = ["Assembly/KmerGenie"]
  end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  end
  def set_default_parameters
    @params['paired'] = dataset_has_column?('Read2')
  end
  def next_dataset
    {'Name'=>'Kmergenie_Result',
     'Report [Link]'=>File.join(@result_dir, "kmergenie_results/histograms_report.html"),
     'Results [File]'=>File.join(@result_dir, "kmergenie_results")
    }
  end
  def commands
    command = "mkdir kmergenie_results\n"
    model = if @params['model'] == 'diploid'
              '--diploid'
            else
              ''
            end
    if  num_reads = @params['sampling'] and num_reads == 'all'
      @dataset_hash.each do |hash|
        hash.each do |k, v|
          if k =~ /\[File/
            command << "echo '#{File.join(@gstore_dir, v)}' >> read_file.txt\n"
          end
        end
      end
    else
      num_lines = num_reads.to_i * 4
      @dataset_hash.each do |hash|
        hash.each do |k, v|
          if k =~ /\[File/
            out_fastq = File.basename(v).gsub(/.gz/, '')
            command << "zcat #{File.join(@gstore_dir, v)}|sort -R|tail -#{num_lines} > #{out_fastq}\n"
            command << "echo '#{out_fastq}' >> read_file.txt\n"
          end
        end
      end
    end
    command << "kmergenie read_file.txt #{model} -t #{@params['cores']} #{@params['cmdOptions']}\n"
    command << "mv read_file.txt kmergenie_results\n"
    command << "mv histograms* kmergenie_results/\n"
  end
end


