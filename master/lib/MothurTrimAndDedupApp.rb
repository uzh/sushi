#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-095604'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class SubreadsApp < SushiFabric::SushiApp
def initialize
super
@name = 'MothurDataClean'
@analysis_category = 'Metagenomics'
@description =<<-EOS
Subreads from SMRT Analysis
<a href='https://mothur.org/wiki/MiSeq_SOP'>https://mothur.org/wiki/MiSeq_SOP</a>
  EOS

@required_columns = ['Name', 'Read1', 'Read2']
@required_params = ['cut-off']
@params['cut-off'] = '0.03'
@params['cut-off', 'description'] = 'Cut-off similarity to cluster OTUss'
@inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
@modules = ["Dev/R"]
end
def next_dataset
report_link_1 = File.join(@result_dir, @dataset['Name'].to_s)
report_link = File.join(report_link_1, 'index.html')
{'Name'=>@dataset['Name'],
#  'Read1 [File]'=>File.join(@result_dir, "#{@dataset['Name']}.filtered_subreads.fastq.gz")
  'Static Report [Link]'=>report_link,
  'SubreadsOut [File]'=>File.join(@result_dir, "#{@dataset['Name']}")
}.merge(extract_columns(@inherit_tags))

end
def commands
run_RApp("EzAppMothurDataClean")
end
end

if __FILE__ == $0
run SubreadsApp
#usecase = Bowtie2App.new

#usecase.project = "p1001"
#usecase.user = 'masamasa'

# set user parameter
# for GUI sushi
#usecase.params['process_mode'].value = 'SAMPLE'
#usecase.params['refBuild'] = 'mm10'
#usecase.params['paired'] = true
#usecase.params['strandMode'] = 'both'
#usecase.params['cores'] = 8
#usecase.params['node'] = 'fgcz-c-048'

# also possible to load a parameterset csv file
# mainly for CUI sushi
#usecase.parameterset_tsv_file = 'tophat_parameterset.tsv'
#usecase.parameterset_tsv_file = 'test.tsv'

# set input dataset
# mainly for CUI sushi
#usecase.dataset_tsv_file = 'tophat_dataset.tsv'

# also possible to load a input dataset from Sushi DB
#usecase.dataset_sushi_id = 3

# run (submit to workflow_manager)
#usecase.run
#usecase.test_run

end
