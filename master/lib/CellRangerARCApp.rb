#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class CellRangerARCApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'CellRangerARCCount'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
This wrapper runs <a href='https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/using/count',>cellranger-arc count</a> in Single-library analysis mode.
     EOS
    @required_columns = ['Name','RNADataDir','ATACDataDir','Species']
    @required_params = ['name', 'refBuild']
    @params['cores'] = ['8', '12', '16']
    @params['cores', "context"] = "slurm"
    @params['ram'] = ['60', '40', '80']
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '300'
    @params['scratch', "context"] = "slurm"
    @params['name'] = 'CellRangerARCCount'
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "referfence genome assembly"
    @params['refFeatureFile'] = 'genes.gtf'
    @params['refFeatureFile', "context"] = "CellRangerARC"
    @params['featureLevel'] = 'gene'
    @params['featureLevel', "context"] = "CellRangerARC"
    @params["excludeIntrons"] = false
    @params['excludeIntrons', 'description'] = 'set to true to reproduce the default behavior in cell ranger v6 and earlier'
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
    @params['keepBam'] = true
    @params['keepBam', 'description'] = 'Keep bam file produced by CellRanger? Usually it is not neccessary for downstream analyses'
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'specify the commandline options for CellRanger (e.g. --include-introns for single nuclei data); do not specify any option that is already covered by the dedicated input fields'
    @params['cmdOptions', "context"] = "CellRangerARC"
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Tools/seqtk", "Dev/R/4.4.0", "Dev/Python", "Tools/samtools"]
    @params['CellRangerARCVersion'] = ["Aligner/CellRangerARC/2.1.0", "Aligner/CellRangerARC/2.0.2", "Aligner/CellRangerARC/2.0.0"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  
  def set_default_parameters
  end
  

  def next_dataset
    report_dir = File.join(@result_dir, "#{@dataset['Name']}_ARC")
    dataset = {
      'Name' => "#{@dataset['Name']}",
      'Species' => @dataset['Species'],
      'refBuild' => @params['refBuild'],
      'refFeatureFile' => @params['refFeatureFile'],
      'featureLevel' => @params['featureLevel'],
      'ResultDir [File]' => report_dir,
      'Report [Link]' => File.join(report_dir, 'web_summary.html'),
      'CountMatrix [Link]' => File.join(report_dir, 'filtered_feature_bc_matrix'),
      'UnfilteredCountMatrix [Link]' => File.join(report_dir, 'raw_feature_bc_matrix'),
      'Read Count' => @dataset['Read Count']
    }.merge(extract_columns(@inherit_tags))
    dataset
  end
  def commands
    command = "module load #{@params["CellRangerARCVersion"]}\n"
    command << run_RApp("EzAppCellRangerARC")
  end
end

if __FILE__ == $0

end
