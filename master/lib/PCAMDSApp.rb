
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class PCAMDSApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'PCAMDS'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Variants'
    @description =<<-EOS
vcf-stats<br/>
    EOS
    @required_columns = ['Name', 'Filtered VCF', 'Grouping File']
    @required_params = ['name']
    @params['cores'] = '1'
    @params['ram'] = '50'
    @params['scratch'] = '100'
    @params['name'] = 'pca_mds'
    @params['mail'] = ""
    @params['Rversion'] = ["Dev/R/4.2.2"]
    @modules = ["Tools/PLINK/1.9beta6.21"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@params['name'],
     'Report [File]'=>report_file,
     'Static Report [Link]'=>report_link,
     'Interactive report [Link]'=>"https://fgcz-shiny.uzh.ch/PopGen_Structure?data=#{report_file}",
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppPCAMDS", lib_path:  "/srv/GT/analysis/jonas/R_LIBS")
    #command = "vcf-stats #{File.join("$GSTORE_DIR", @dataset[0]['Filtered VCF [File]'])} -p #{@params['name']}/vcf_stats"
  end
end
