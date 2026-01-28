#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class SpaceRangerApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'SpaceRangerCount'
    @analysis_category = 'Spatial'
    @description =<<-EOS
This wrapper runs <a href='https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/using/count',>space ranger count</a> in Single-library analysis mode.
    EOS
    @required_columns = ['Name','RawDataDir','Species','Slide','Area']
    @required_params = ['name', 'refBuild']
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '60'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '200'
    @params['scratch', "context"] = "slurm"
    @params['name'] = 'SpaceRangerCount'
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "referfence genome assembly"
    @params['refFeatureFile'] = 'genes.gtf'
    @params['refFeatureFile', "context"] = "SpaceRanger"
    @params['featureLevel'] = 'gene'
    @params['featureLevel', "context"] = "SpaceRanger"
    @params['transcriptTypes'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA', 'long_noncoding', 'short_noncoding', 'pseudogene']
    @params['transcriptTypes', 'multi_selection'] = true
    @params['transcriptTypes', 'selected'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA']
    @params['probesetFile'] =  {'select'=>''}
    Dir["/srv/GT/databases/10x_Probesets/Visium/*"].sort.select{|design| File.file?(design)}.each do |dir|
      @params['probesetFile'][File.basename(dir)] = File.basename(dir)
    end
    @params['customProbesFile'] = ''
    @params['customProbesFile', 'file_upload'] = true
    @params['customProbesFile', 'description'] = 'Custom probeset CSV-file according to 10x specifications (https://tinyurl.com/VisiumProbeSetsDef).Note that all genes listed must have a corresponding entry in secondRef or controlSeqs. Custom probes must have the same length as the probes in the reference file.'

    @params['panelFile'] =  {'select'=>''}
    Dir["/srv/GT/databases/10x/Visium/panels/*"].sort.select{|design| File.file?(design)}.each do |dir|
      @params['panelFile'][File.basename(dir)] = File.basename(dir)
    end
    @params['panelFile', 'description'] = 'for protein panels'
    @params['controlSeqs'] = ''
    @params['controlSeqs', 'description'] = 'The extra control sequences (such as spikein sequences) available in https://fgcz-gstore.uzh.ch/reference/controlSeqs.fa'
    @params['secondRef'] = ''
    @params['secondRef', 'description'] = 'full path to fasta file with e.g. viralGenes'
    @params['keepAlignment'] = true
    @params['keepAlignment', 'description'] = 'Keep bam/cram file produced by SpaceRanger? Usually it is not neccessary for downstream analyses'
    @params['runSegmentation'] = false
    @params['runSegmentation', 'description'] = '10x built in segmentation for H&E images, currently very experimental and only available for >=v4'
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'specify the commandline options for SpaceRanger; do not specify any option that is already covered by the dedicated input fields'
    @params['cmdOptions', "context"] = "SpaceRanger"
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @params['SpaceRangerVersion'] = ["Aligner/SpaceRanger/4.0.1","Aligner/SpaceRanger/3.1.3","Aligner/SpaceRanger/3.1.2","Aligner/SpaceRanger/3.0.1"]
    @modules = ["Dev/R", "Aligner/CellRanger", "Tools/samtools"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
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
        'CountMatrix [Link]'=>File.join(report_dir, 'filtered_feature_bc_matrix'),
        'Read Count'=>@dataset['Read Count'],
        'Anndata [Link]' => (
          @params['SpaceRangerVersion'] == "Aligner/SpaceRanger/4.0.1" ?
            File.join(report_dir, 'segmented_outputs/filtered_feature_cell_matrix.h5') :
            File.join(report_dir, 'binned_outputs/square_008um/filtered_feature_bc_matrix.h5')
        ),
        'SourceImage [Link]'=>@dataset['Image'],
        'BinnedOutput [Link]'=>File.join(report_dir, 'binned_outputs/square_002um'),
        'SpaceRangerDir [Link]'=>report_dir,
        'Count [Link]'=>File.join(report_dir, "#{@dataset['Name']}-counts.txt")        
      }.merge(extract_columns(@inherit_tags))
   if @params['keepAlignment']
       dataset['AlignmentFile [Link]'] = File.join(report_dir, 'possorted_genome_bam.cram')
    end    
    dataset
  end
  def commands
    command = "module load  #{@params["SpaceRangerVersion"]}\n"
    command << run_RApp("EzAppSpaceRanger")
  end
end

if __FILE__ == $0

end
