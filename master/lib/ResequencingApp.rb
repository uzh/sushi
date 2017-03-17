#!/usr/bin/env ruby
# encoding: utf-8
Version = '20170313'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class ResequencingApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Resequencing'
    @analysis_category = 'Map'
    @description =<<-EOS
Resequencing pipeline (mapping, variant calling and consensusing from SMRT Analysis
<a href='https://github.com/PacificBiosciences/SMRT-Analysis/wiki'>https://github.com/PacificBiosciences/SMRT-Analysis/wiki</a>
EOS
    
    @required_columns = ['Name','Reads']
    @required_params = ['RefFile']
    @params['RefFile'] = {'select'=>''}
    @params['RefFile'][lambda] = "common/userdata.d/references/lambda"
    @params['RefFile'][ecoli] = "common/userdata.d/references/ecoli"
    @params['RefFile'][bsubtilis] = "common/userdata.d/references/B_subtilis_v1"
    @params['RefFile'][yeast] = "common/userdata.d/references/s_cerevisiae"
    @params['RefFile'][yhalotolerans] = "common/userdata.d/references/yhalotolerans"
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '50'
    @params['scratch'] = '400'
    @params['minSubReadLength'] = '50'
    @params['minSubReadLength', 'description'] = 'Minimum Subread Length. Subreads shorter than this value (in bps) are filtered out and excluded from pre-assembly.'
    @params['readScore'] = '75'
    @params['readScore', 'description'] = 'Minimum Polymerase Read Quality. Polymerase reads with lower quality than this value are filtered out and excluded from pre-assembly.'
    @params['minLength'] = '50'
    @params['minLength', 'description'] = 'Minimum Polymerase Read Length. Polymerase reads shorter than this value (in bps) are filtered out and excluded from pre-assembly.'
    @params['maxDivergence'] = '30'
    @params['maxDivergence', 'description'] = 'The maximum allowed divergence (in %) of a read from the reference sequence.'
    @params['minAnchorSize'] = '12'
    @params['minAnchorSize', 'description'] = 'The minimum size of the read (in bps) that must match against the reference.'
    @params['mail'] = ""
  end
  def next_dataset
    {'Name'=>@dataset['Name'],  
     'BAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam"), 
     'BAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam.bai"),
     'VCF [File]'=>File.join(@result_dir, "#{@dataset['Name']}.vcf.gz"),
     'VCFINDEX [File]'=>File.join(@result_dir, "#{@dataset['Name']}.vcf.gz.tbi"),
     'Consensus [File]'=>File.join(@result_dir, "#{@dataset['Name']}.polished.fasta"),	
     'ResequencingOut [File]'=>File.join(@result_dir, "#{@dataset['Name']}"),  
  }.merge(extract_column("Factor")).merge(extract_column("B-Fabric"))
  end
  def commands
    run_RApp("EzAppResequencing")
  end
end

if __FILE__ == $0
  run ResequencingApp
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

