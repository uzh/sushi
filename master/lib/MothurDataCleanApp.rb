#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-095604'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class MothurDataCleanApp < SushiFabric::SushiApp
def initialize
super
@name = 'MothurDataClean'
@analysis_category = 'Metagenomics'
@description =<<-EOS
OTU-based metagenomics analysis with Mothur.
<a href='https://mothur.org/wiki/MiSeq_SOP'>https://mothur.org/wiki/MiSeq_SOP</a>
  EOS
@params['process_mode'] = 'DATASET'
@required_columns = ['Name', 'Read1', 'Technology']
@required_params = ['cutOff', 'diffs_Illumina','diffs_PacBio']
@params['cores'] = '1'
@params['ram'] = '8'
@params['scratch'] = '10'
@params['cutOff'] = '80'
@params['cutOff', 'description'] = 'Cut-off for taxonomy assignment in classify.seqs'
@params['diffs_Illumina'] = '2'
@params['diffs_Illumina', 'description'] = 'Differences allowed in the pre.cluster step. Should be 1 every 100 bases.If the data are only Pacbio, it is ignored'
@params['diffs_PacBio'] = '15'
@params['diffs_PacBio', 'description'] = 'Differences allowed in the pre.cluster step. Should be 1 every 100 bases.If the data are only Illumina, it is ignored'
@params['minLen_Illumina'] = '145'
@params['minLen_Illumina', 'description'] = 'Illumina sequences shorter than this long are removed.'
@params['maxLen_Illumina'] = '155'
@params['maxLen_Illumina', 'description'] = 'Illumina sequences longer than this are removed.'
@params['minLen_PacBio'] = '1400'
@params['minLen_PacBio', 'description'] = 'PacBio sequences shorter than this long are removed.'
@params['maxLen_PacBio'] = '1700'
@params['maxLen_PacBio', 'description'] = 'PacBio sequences longer than this are removed.'
@params['mail'] = ""
@inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
@modules = ["Dev/R"]
end
def next_dataset
@params['name'] = "MothurDataClean"
report_file = File.join(@result_dir, '00index_files')
report_link = File.join(@result_dir, '00index.html')
{'Name'=>@params['name'],
  'CountTablePacbio [File]'=>File.join(@result_dir, "PacBio.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table"),
  'PreClusteredFastaFilePacbio [File]'=>File.join(@result_dir, "PacBio.good.unique.good.filter.unique.precluster.pick.pick.fasta"),
  'TaxonomyFilePacbio [File]'=>File.join(@result_dir, "PacBio.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy"),
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
