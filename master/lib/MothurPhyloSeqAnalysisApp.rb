#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-095604'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class MothurPhyloSeqAnalysisApp < SushiFabric::SushiApp
def initialize
super
@name = 'MothurPhyloSeqAnalysis'
@analysis_category = 'Metagenomics'
@description =<<-EOS
16S metagenmics visualization with Phyloseq.
<a href='http://joey711.github.io/phyloseq/index.html'>http://joey711.github.io/phyloseq/index.html</a>
  EOS
@params['process_mode'] = 'DATASET'
@required_columns = ['Name', 'OTU_pacbio','Taxonomy_pacbio','OTU_Illumina','Taxonomy_Illumina']
@required_params = ['designMatrix']
@params['cores'] = '1'
@params['ram'] = '8'
@params['scratch'] = '10'
@params['designMatrix'] = ''
@params['designMatrix', 'description'] = 'Tab-separated file assigning samples to groups.'
@params['mail'] = ""
@inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
@modules = ["Dev/R"]
end
def next_dataset
@params['name'] = "MothurPhyloSeq"
report_file = File.join(@result_dir, '00index_files')
report_link = File.join(@result_dir, '00index.html')
{'Name'=>@dataset['Name'],
  'Report [File]'=>report_file,
  'Static Report [Link]'=>report_link
}

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
