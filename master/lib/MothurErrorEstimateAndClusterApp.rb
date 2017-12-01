#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-095604'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class MothurErrorEstimateAndClusterApp < SushiFabric::SushiApp
def initialize
super
@name = 'MothurErrorEstimateAndCluster'
@analysis_category = 'Metagenomics'
@description =<<-EOS
OTU-based metagenomics analysis with Mothur.
<a href='https://mothur.org/wiki/MiSeq_SOP'>https://mothur.org/wiki/MiSeq_SOP</a>
  EOS
@params['process_mode'] = 'DATASET'
@required_columns = ['Name', 'CountTablePacbio','PreClusteredFastaFilePacbio','TaxonomyFilePacbio', 'CountTableIllumina','PreClusteredFastaFileIllumina', 'TaxonomyFileIllumina',]
@required_params = ['cutOff','referenceGroup', 'referenceFasta']
@params['cores'] = '1'
@params['ram'] = '8'
@params['scratch'] = '10'
@params['cutOff'] = '0.03'
@params['cutOff', 'description'] = 'Cut-off similarity to cluster OTUs'
@params['referenceGroup'] = 'Mock'
@params['referenceGroup', 'description'] = 'Group of sample against which to estimate error rates'
@params['referenceFasta'] = ''
@params['referenceFasta', 'description'] = 'Full path to fasta file for the mock community.'
@params['mail'] = ""
@inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
@modules = ["Dev/R"]
end
def next_dataset
@params['name'] = "MothurErrorEstAndClustering"
report_file = File.join(@result_dir, '00index_files')
report_link = File.join(@result_dir, '00index.html')
pacbioTaxFile = "PacBio.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc." + @params['cutOff'] + ".cons.taxonomy"
illuminaTaxFile = "Illumina.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc." + @params['cutOff'] + ".cons.taxonomy"
{'Name'=>@dataset['Name'],
  'OTU_pacbio [File]'=>File.join(@result_dir, "PacBio.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared"),
  'Taxonomy_pacbio [File]'=>File.join(@result_dir, pacbioTaxFile),
  'OTU_Illumina [File]'=>File.join(@result_dir, "Illumina.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared"),
  'Taxonomy_Illumina [File]'=>File.join(@result_dir, illuminaTaxFile),
  'Report [File]'=>report_file,
  'Static Report [Link,File]'=>report_link,
}

end
def commands
run_RApp("EzAppMothurErrorEstimateAndCluster")
end
end

if __FILE__ == $0

end
