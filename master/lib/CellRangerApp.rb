#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class CellRangerApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'CellRangerCount'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
This wrapper runs <a href='https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count',>cellranger count</a> in Single-library analysis mode.
     EOS
    @required_columns = [['Name','RawDataDir','Species'], ['Name', 'Read1', 'Read2', 'Species']]
    @required_params = ['name', 'refBuild']
    @params['cores'] = ['8', '12', '16']
    @params['cores', "context"] = "slurm"
    @params['ram'] = ['60', '40', '80']
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '300'
    @params['scratch', "context"] = "slurm"
    @params['name'] = 'CellRangerCount'
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "referfence genome assembly"
    @params['refFeatureFile'] = 'genes.gtf'
    @params['refFeatureFile', "context"] = "CellRanger"
    @params['featureLevel'] = 'gene'
    @params['featureLevel', "context"] = "CellRanger"
    @params['TenXLibrary'] = ['GEX', 'VDJ']
    @params['TenXLibrary', 'description'] = 'Which 10X library?'
    @params['chemistry'] = ['auto', 'threeprime', 'fiveprime', 'SC3Pv1', 'SC3Pv2', 'SC3Pv3', 'SC5P-PE', 'SC5P-R2', 'ARC-v1']
    @params['chemistry', 'description'] = 'Assay configuration. NOTE: by default the assay configuration is detected automatically, which is the recommended mode.'
    @params["includeIntrons"] = true
    @params['includeIntrons', 'description'] = 'set to false to reproduce the default behavior in cell ranger v6 and earlier'
    @params['expectedCells'] = ''
    @params['expectedCells', 'description'] = 'Expected number of recovered cells. Leave this free to let cellranger estimate the cell number.'
    @params['transcriptTypes'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA', 'long_noncoding', 'short_noncoding', 'pseudogene']
    @params['transcriptTypes', 'multi_selection'] = true
    @params['transcriptTypes', 'selected'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA']
    @params['controlSeqs'] = ''
    @params['controlSeqs', 'description'] = 'The extra control sequences (such as spikein sequences) available in https://fgcz-gstore.uzh.ch/reference/controlSeqs.fa'
    @params['secondRef'] = ''
    @params['secondRef', 'description'] = 'full path to fasta file with e.g. viralGenes'
    @params['runVeloCyto'] = false
    @params['runVeloCyto', 'description'] = 'generate loom file by velocyto?'
    @params['bamStats'] = false
    @params['bamStats', 'description'] = 'Compute stats per cell from the bam file?'
    @params['keepAlignment'] = true
    @params['keepAlignment', 'description'] = 'Keep cram/bam file produced by CellRanger? Usually it is not neccessary for downstream analyses'
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'specify the commandline options for CellRanger (e.g. --include-introns for single nuclei data); do not specify any option that is already covered by the dedicated input fields'
    @params['cmdOptions', "context"] = "CellRanger"
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Tools/seqtk", "Dev/R", "Dev/Python/3.11.5_dev", "Tools/samtools"]
    @params['CellRangerVersion'] = ["Aligner/CellRanger/10.0.0", "Aligner/CellRanger/9.0.0", "Aligner/CellRanger/8.0.1", "Aligner/CellRanger/7.1.0"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def set_default_parameters
  end
  def next_dataset
    report_dir = File.join(@result_dir,"#{@dataset['Name']}")
    if @params['TenXLibrary'] == "VDJ"
      dataset = {
        'Name'=>@dataset['Name'],
        'Species'=>@dataset['Species'],
        'refBuild'=>@params['refBuild'],
        'refFeatureFile'=>@params['refFeatureFile'],
        'featureLevel'=>@params['featureLevel'],
        'ResultDir [File]'=>report_dir,
        'Report [Link]'=>File.join(report_dir, 'web_summary.html'),
        'Read Count'=>@dataset['Read Count']
      }
    else
      dataset = {
        'Name'=>@dataset['Name'],
        'Species'=>@dataset['Species'],
        'refBuild'=>@params['refBuild'],
        'refFeatureFile'=>@params['refFeatureFile'],
        'featureLevel'=>@params['featureLevel'],
        'transcriptTypes'=>@params['transcriptTypes'],
        'SCDataOrigin'=>'10X',
        'ResultDir [File]'=>report_dir,
        'Report [Link]'=>File.join(report_dir, 'web_summary.html'),
        'CountMatrix [Link]'=>File.join(report_dir, 'filtered_feature_bc_matrix'),
        'UnfilteredCountMatrix [Link]'=>File.join(report_dir, 'raw_feature_bc_matrix'),
        'Read Count'=>@dataset['Read Count']
      }
      if @params['keepAlignment']
       dataset['AlignmentFile [Link]'] = File.join(report_dir, 'possorted_genome_bam.cram')
      end
    end
    dataset.merge(extract_columns(@inherit_tags))
  end
  def commands
    command = "module load  #{@params["CellRangerVersion"]}\n"
    command << run_RApp("EzAppCellRanger")
  end
end

if __FILE__ == $0

end
