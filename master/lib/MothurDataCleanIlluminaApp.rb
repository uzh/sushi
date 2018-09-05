#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-095604'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class MothurDataCleanIlluminaApp < SushiFabric::SushiApp
def initialize
super
@name = 'MothurDataCleanIllumina'
@analysis_category = 'Metagenomics'
@description =<<-EOS
OTU-based metagenomics analysis of Illumina data with Mothur.
<a href='https://mothur.org/wiki/MiSeq_SOP'>https://mothur.org/wiki/MiSeq_SOP</a>
  EOS
@params['process_mode'] = 'DATASET'
@required_columns = ['Name', 'Read1']
@required_params = ['cutOff', 'diffs_Illumina']
@params['cores'] = '1'
@params['ram'] = '8'
@params['scratch'] = '10'
@params['cutOff'] = '80'
@params['cutOff', 'description'] = 'Cut-off for taxonomy assignment in classify.seqs'
@params['diffs_Illumina'] = '2'
@params['diffs_Illumina', 'description'] = 'Differences allowed in the pre.cluster step. Should be 1 every 100 bases.If the data are only Pacbio, it is ignored'
@params['minLen_Illumina'] = '145'
@params['minLen_Illumina', 'description'] = 'Sequences shorter than this long are removed.'
@params['maxLen_Illumina'] = '330'
@params['maxLen_Illumina', 'description'] = 'Sequences longer than this are removed.'
@params['mail'] = ""
@inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
@modules = ["Dev/R"]
end
def next_dataset
@params['name'] = "MothurDataCleanIllumina"
report_file = File.join(@result_dir, '00index_files')
report_link = File.join(@result_dir, '00index.html')
{'Name'=>@params['name'],
  'CountTableIllumina [File]'=>File.join(@result_dir, "Illumina.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table"),
  'PreClusteredFastaFileIllumina [File]'=>File.join(@result_dir, "Illumina.good.unique.good.filter.unique.precluster.pick.pick.fasta"),
  'TaxonomyFileIllumina [File]'=>File.join(@result_dir, "Illumina.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy"),
  'Report [File]'=>report_file,
  'Static Report [Link,File]'=>report_link,
  }

end
def commands
run_RApp("EzAppMothurDataClean")
end
end

if __FILE__ == $0
end
