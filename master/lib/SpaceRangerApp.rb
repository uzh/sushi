#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class SpaceRangerApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'SpaceRangerCount'
    @analysis_category = 'SpatialTrx'
    @description =<<-EOS
This wrapper runs <a href='https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/using/count',>space ranger count</a> in Single-library analysis mode.
    EOS
    @required_columns = ['Name','RawDataDir','Species', 'Image', 'Slide', 'Area']
    @required_params = ['name', 'refBuild']
    @params['cores'] = '8'
    @params['ram'] = '60'
    @params['scratch'] = '200'
    @params['name'] = 'SpaceRangerCount'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['featureLevel'] = 'gene'
    @params['transcriptTypes'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA', 'long_noncoding', 'short_noncoding', 'pseudogene']
    @params['transcriptTypes', 'multi_selection'] = true
    @params['transcriptTypes', 'selected'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA']
    @params['controlSeqs'] = ''
    @params['controlSeqs', 'description'] = 'The extra control sequences (such as spikein sequences) available in https://fgcz-gstore.uzh.ch/reference/controlSeqs.fa'
    @params['probesetFile'] =  {'select'=>''}
    Dir["/srv/GT/databases/10x_Probesets/Visium/*"].sort.select{|design| File.file?(design)}.each do |dir|
      @params['probesetFile'][File.basename(dir)] = File.basename(dir)
    end
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'specify the commandline options for SpaceRanger; do not specify any option that is already covered by the dedicated input fields'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @params['SpaceRangerVersion'] = ["Aligner/SpaceRanger/2.1.0", "Aligner/SpaceRanger/2.0.1"]
    @modules = ["Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def set_default_parameters
  end
  def next_dataset
    report_dir = File.join(@result_dir,"#{@dataset['Name']}")
    dataset = {
        'Name'=>@dataset['Name'],
        'Species'=>@dataset['Species'],
        'refBuild'=>@params['refBuild'],
        'refFeatureFile'=>@params['refFeatureFile'],
        'featureLevel'=>@params['featureLevel'],
        'transcriptTypes'=>@params['transcriptTypes'],
        'ResultDir [File]'=>report_dir,
        'Report [Link]'=>File.join(report_dir, 'web_summary.html'),
        'CountMatrix [Link]'=>File.join(report_dir, 'filtered_feature_bc_matrix')
      }.merge(extract_columns(@inherit_tags))
    dataset
  end
  def commands
    command = "module load  #{@params["SpaceRangerVersion"]}\n"
    command << run_RApp("EzAppSpaceRanger")
  end
end

if __FILE__ == $0

end
