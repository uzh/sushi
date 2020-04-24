#!/usr/bin/env ruby
# encoding: utf-8
Version = '20200424-195259'

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
    @inherit_tags = ["Species", "refBuild"]
    @modules = ["Variants/Beagle/5.1"]
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'Phased and Imputed VCF [File]'=>File.join(@result_dir, "#{@dataset['Name']}_phased_imputed.vcf.gz"),
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    filtered_vcf = File.join(@gstore_dir, @dataset['Filtered VCF'])
    command = "java -Xmx#{@params['ram']}g -jar $Beagle_jar gt=#{filtered_vcf} out=#{@dataset['Name']}_phased_imputed nthreads=#{@params['cores']}\n"
  end
end


