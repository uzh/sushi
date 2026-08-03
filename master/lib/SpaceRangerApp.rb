#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class SpaceRangerApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'SpaceRangerCount'
    # spaceranger count unconditional. samtools (CRAM) gated on keepAlignment.
    # Biostrings/rtracklayer gated on controlSeqs (via shared
    # getCellRangerGEXReference helper) -- this app does not declare secondRef
    # or extendThreePrime, so only the controlSeqs-gated path is reachable.
    @citation = [
      '10x Genomics. Space Ranger. https://www.10xgenomics.com/support/software/space-ranger',
      'Ståhl, P.L. et al. Visualization and analysis of gene expression in tissue sections by spatial transcriptomics. Science 353(6294), 78-82 (2016). https://doi.org/10.1126/science.aaf2403',
      'Li, H. et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics 25(16), 2078-2079 (2009). https://doi.org/10.1093/bioinformatics/btp352',
      'Pagès, H., Aboyoun, P., Gentleman, R. & DebRoy, S. Biostrings: Efficient manipulation of biological strings. R package version 2.80.1. https://doi.org/10.18129/B9.bioc.Biostrings',
      'Lawrence, M., Gentleman, R. & Carey, V. rtracklayer: an R package for interfacing with genome browsers. Bioinformatics 25(14), 1841-1842 (2009). https://doi.org/10.1093/bioinformatics/btp328'
    ]
    @analysis_category = 'Spatial'
    @description =<<-EOS
This wrapper runs <a href='https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/using/count',>space ranger count</a> in Single-library analysis mode.
    EOS
    @required_columns = [['Name','RawDataDir','Species','Slide','Area'], ['Name','Read1','Read2','Species','Slide','Area']]
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
    @params['probesetFile', 'description'] = 'Required for probe-based assays (FFPE, CytAssist, VisiumHD). Not needed for standard Visium frozen tissue (3\' polyA capture). Probe set compatibility: v2.x = VisiumHD, CytAssist (FFPE/FF), CytAssist Gene+Protein; v1.x = original Visium FFPE (SD), Mouse CytAssist (FFPE/FF/Fixed Frozen). Recommend v2.1.0 for latest reference genomes (GRCh38-2024-A / GRCm39-2024-A).'
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
    @params['keepAlignment', 'description'] = 'Keep bam/cram file produced by SpaceRanger? Usually not necessary for downstream analyses.'
    @params['runSegmentation'] = false
    @params['runSegmentation', 'description'] = '10x built-in nuclear segmentation for H&E images. Applies to VisiumHD runs only (H-slides). Currently experimental, requires SpaceRanger >= v4.'
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'Specify commandline options for SpaceRanger; do not specify options already covered by dedicated input fields. Use --unknown-slide for VisiumHD slides not in the 10x database.'
    @params['cmdOptions', "context"] = "SpaceRanger"
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @params['SpaceRangerVersion'] = ["Aligner/SpaceRanger/4.1.0","Aligner/SpaceRanger/4.0.1","Aligner/SpaceRanger/3.1.3","Aligner/SpaceRanger/3.1.2","Aligner/SpaceRanger/3.0.1"]
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
        'CountMatrix [Link]'=>File.join(report_dir, 'filtered_feature_bc_matrix'),
        'Read Count'=>@dataset['Read Count'],
        'SourceImage [Link]'=>@dataset['Image'],
        'SpaceRangerDir [Link]'=>report_dir,
        'Count [Link]'=>File.join(report_dir, "#{@dataset['Name']}-counts.txt")
      }.merge(extract_columns(@inherit_tags))
    if @dataset['Slide'] && @dataset['Slide'].start_with?('H')
      dataset['BinnedOutput [Link]'] = File.join(report_dir, 'binned_outputs/square_002um')
    end
    if @params['keepAlignment']
      dataset['AlignmentFile [Link]'] = File.join(report_dir, 'possorted_genome_bam.cram')
    end
    if @params['runSegmentation']
      dataset['Anndata [Link]'] = File.join(report_dir, 'segmented_outputs/filtered_feature_cell_matrix.h5')
    end
    # SUSHI renders columns in REVERSE hash-insertion order (after Name & Factor sort).
    # Insert Report last so it appears as the leftmost link column in the UI.
    dataset['Report [Link]'] = File.join(report_dir, 'web_summary.html')
    dataset
  end
  def commands
    command = "module load  #{@params["SpaceRangerVersion"]}\n"
    command << run_RApp("EzAppSpaceRanger")
  end
end

if __FILE__ == $0

end
