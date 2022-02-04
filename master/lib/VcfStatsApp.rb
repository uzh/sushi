#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class VcfStatsApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'VcfStats'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Variants'
    @description =<<-EOS
vcf-stats<br/>
    EOS
    @required_columns = ['Name', 'Filtered VCF', 'Dummy']
    @required_params = ['name']
    @params['cores'] = '1'
    @params['ram'] = '50'
    @params['scratch'] = '100'
    @params['name'] = 'vcf_stats'
    @params['mail'] = ""
    @modules = ["Tools/vcftools/0.1.16"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@params['name'],
     'Report [File]'=>report_file,
     'Html [Link]'=>report_link,
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppVcfStats")
    #command = "vcf-stats #{File.join("$GSTORE_DIR", @dataset[0]['Filtered VCF [File]'])} -p #{@params['name']}/vcf_stats"
  end
end
