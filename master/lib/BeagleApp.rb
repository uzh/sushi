#!/usr/bin/env ruby
# encoding: utf-8
Version = '20200424-200202'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class BeagleApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Beagle'
    @description =<<-EOS
Beagle is a software package for phasing genotypes and for imputing ungenotyped markers.<br />
<br />
* https://faculty.washington.edu/browning/beagle/beagle.html 
    EOS

    @analysis_category = 'Polyploid'
    @required_columns = ['Name','Filtered VCF']
    @required_params = []
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '50'
    @params['scratch'] = '100'
    @params['DataSet'] = []
    @inherit_tags = ["Species", "refBuild"]
    @modules = ["Variants/Beagle/5.1"]
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'Phased and Imputed VCF [File]'=>File.join(@result_dir, "#{@dataset['Name']}_phased_imputed.vcf.gz"),
    }.merge(extract_columns(@inherit_tags))
  end
  def set_default_parameters
    if data_set = DataSet.find_by_id(@dataset_sushi_id)
      @params['DataSet'] = Hash[*data_set.project.data_sets.map{|d| [d.id.to_s + ":" + d.name, d.id]}.sort_by{|name, id| id}.reverse.flatten]
    end
  end
  def sample_path(data_set)
    paths = []
    data_set.samples.each do |sample|
      sample.to_hash.each do |header, file|
        if (header.tag?('File') or header.tag?('Link')) and file !~ /http/
          paths << File.dirname(file)
        end
      end
    end
    paths.uniq!
    if path = paths.first
      path.split('/')[0,2].join('/')
    end
  end
  def commands
    filtered_vcf = File.join(@gstore_dir, @dataset['Filtered VCF'])
    if ref_dataset = DataSet.find_by_id(params['DataSet']) and
       ref_vcf_gz_dir = File.join(@gstore_dir, sample_path(ref_dataset)) and
       ref_vcf_gz = Dir[File.join(ref_vcf_gz_dir, "*.vcf.gz")].to_a.first

      command = "java -Xmx#{@params['ram']}g -jar $Beagle_jar ref=#{ref_vcf_gz} gt=#{filtered_vcf} out=#{@dataset['Name']}_phased_imputed nthreads=#{@params['cores']}\n"
    else
      command = "java -Xmx#{@params['ram']}g -jar $Beagle_jar gt=#{filtered_vcf} out=#{@dataset['Name']}_phased_imputed nthreads=#{@params['cores']}\n"
    end

  end
end


