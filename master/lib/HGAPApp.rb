#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-094629'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class HGAPApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'HGAP'
    @analysis_category = 'Assemble'
    @description =<<-EOS
The Hierarchical Genome Assembly Process (HGAP) from SMRT Analysis
<a href='https://github.com/PacificBiosciences/SMRT-Analysis/wiki'>https://github.com/PacificBiosciences/SMRT-Analysis/wiki</a>
EOS
    
    @required_columns = ['Name','Reads']
    @required_params = ['genomeSize', 'xCoverage']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '300'
    @params['minSubReadLength'] = '500'
    @params['minSubReadLength', 'description'] = 'Minimum Subread Length. Subreads shorter than this value (in bps) are filtered out and excluded from pre-assembly.'
    @params['readScore'] = '0.8'
    @params['readScore', 'description'] = 'Minimum Polymerase Read Quality. Polymerase reads with lower quality than this value are filtered out and excluded from pre-assembly.'
    @params['minLength'] = '100'
    @params['minLength', 'description'] = 'Minimum Polymerase Read Length. Polymerase reads shorter than this value (in bps) are filtered out and excluded from pre-assembly.'
    @params['genomeSize'] = '5000000'
    @params['genomeSize', 'description'] = 'The approximate genome size, in bps.' 
    @params['xCoverage'] = '25'
    @params['xCoverage', 'description'] = 'Fold coverage to target for when picking the minimum fragment length for assembly; typically 15 to 25.'
    @params['targetChunks'] = '6'
    @params['targetChunks', 'description'] = 'The number of pieces to split the data files into while running pre-assembly.'
    @params['splitBestn'] = '10'
    @params['splitBestn', 'description'] = 'The number of alignments to consider for each read for a particular chunk.'
    @params['totalBestn'] = '24'
    @params['totalBestn', 'description'] = 'The number of potential alignments to consider across all chunks for a particular read.'
    @params['minCorCov'] = '6'
    @params['minCorCov', 'description'] = 'The minimum coverage to maintain correction for a read.  If the coverage falls below that threshold, the read will be broken at that juntion.'
    @params['ovlErrorRate'] = '0.06'
    @params['ovlErrorRate', 'description'] = 'Trimming and assembly overlaps above this error limit are not computed.'
    @params['ovlMinLen'] = '40'
    @params['ovlMinLen', 'description'] = 'Overlaps shorter than this length (in bps) are not computed.'
    @params['merSize'] = '14'
    @params['merSize', 'description'] = 'The length of the seeds (in base pairs) used by the seed-and-extend algorithm.'
    @params['maxDivergence'] = '30'
    @params['maxDivergence', 'description'] = 'The maximum allowed divergence (in %) of a read from the reference sequence.'
    @params['minAnchorSize'] = '12'
    @params['minAnchorSize', 'description'] = 'The minimum size of the read (in bps) that must match against the reference.'
    @params['mail'] = ""
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
    @modules = ["Dev/R"]
  end
  def next_dataset
        report_link_1 = File.join(@result_dir, @dataset['Name'].to_s)
    report_link = File.join(report_link_1, 'index.html')
    {'Name'=>@dataset['Name'], 
     'Reads'=>@dataset['Reads'], 
      'Static Report [Link]'=>report_link,
     'Draft [File]'=>File.join(@result_dir, "#{@dataset['Name']}.contigs.fasta"),	
     'HGAPOut [File]'=>File.join(@result_dir, "#{@dataset['Name']}"),
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppHGAP")
  end
end

if __FILE__ == $0
  run HGAPApp
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

