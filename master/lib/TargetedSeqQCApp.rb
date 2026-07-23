#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class TargetedSeqQCApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'TargetedSeqQC'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'QC'
    @description =<<-EOS
Cohort QC for targeted / exome sequencing (WES).<br/>
Runs Picard CollectHsMetrics, mosdepth, samtools (stats/flagstat/idxstats) and
optionally Picard CollectInsertSizeMetrics / VerifyBamID2 for every aligned BAM,
then aggregates everything into a single <b>MultiQC</b> report.<br/>
A modern, MultiQC-based alternative to the TEQC-based Teqc app.
    EOS
    @required_columns = ['Name', 'BAM', 'BAI', 'refBuild', 'Species']
    @required_params = ['name', 'refBuild', 'designName']
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '50'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '150'
    @params['scratch', "context"] = "slurm"
    @params['name'] = 'TargetedSeqQC'
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "reference genome assembly"
    @params['designName'] = ''
    @params['designName', 'description'] = "capture design: a sub-folder name under #{TEQC_DESIGN_DIR} (e.g. Twist_Bioscience_for_Illumina_Exome_2_5_hg38) or a full path to a target BED"
    @params['designName', "context"] = "TargetedSeqQC"
    @params['vendor'] = ['auto', 'illumina', 'agilent']
    @params['vendor', 'description'] = "illumina: one BED as bait+target; agilent: *_Covered.bed (bait) + *_Regions.bed (target)"
    @params['vendor', "context"] = "TargetedSeqQC"
    @params['restrictToMainChromosomes'] = true
    @params['restrictToMainChromosomes', 'description'] = "restrict QC to assembled chromosomes, dropping scaffolds/decoys"
    @params['restrictToMainChromosomes', "context"] = "TargetedSeqQC"
    @params['mainChromMaxNchar'] = '6'
    @params['mainChromMaxNchar', 'description'] = "a contig is a 'main chromosome' if its name is at most this many characters (chr1..chr22,chrX,chrY,chrM); scaffolds like GL000009.2 are longer"
    @params['mainChromMaxNchar', "context"] = "TargetedSeqQC"
    @params['mapq'] = '20'
    @params['mapq', "context"] = "TargetedSeqQC"
    @params['baseq'] = '20'
    @params['baseq', "context"] = "TargetedSeqQC"
    @params['coverageThresholds'] = '1,10,20,30,50,100'
    @params['coverageThresholds', "context"] = "TargetedSeqQC"
    @params['runInsertSize'] = true
    @params['runInsertSize', "context"] = "TargetedSeqQC"
    @params['markDuplicates'] = false
    @params['markDuplicates', 'description'] = "run Picard MarkDuplicates first (only if the BAMs are NOT already duplicate-marked)"
    @params['markDuplicates', "context"] = "TargetedSeqQC"
    @params['verifyBamIDSVD'] = ''
    @params['verifyBamIDSVD', 'description'] = "path to a VerifyBamID2 SVD resource prefix to enable contamination (FREEMIX); empty to skip"
    @params['verifyBamIDSVD', "context"] = "TargetedSeqQC"
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Tools/samtools", "Tools/mosdepth", "Tools/Picard", "QC/MultiQC", "Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end

  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
  end

  def next_dataset
    report_dir = File.join(@result_dir, @params['name'])
    {
      'Name' => @params['name'],
      'Report [File]' => report_dir,
      'Html [Link]' => File.join(report_dir, '00index.html'),
      'MultiQC Report [Link]' => File.join(report_dir, 'multiqc', 'multiqc_report.html'),
      'Species' => (dataset = @dataset.first and dataset['Species']),
      'refBuild' => @params['refBuild']
    }.merge(extract_columns(@inherit_tags))
  end

  def commands
    run_RApp("EzAppTargetedSeqQC")
  end
end

if __FILE__ == $0
  usecase = TargetedSeqQCApp.new
  usecase.project = "p1001"
  usecase.user = "masa"
  usecase.parameterset_tsv_file = 'tophat_parameterset.tsv'
  usecase.dataset_tsv_file = 'tophat_dataset.tsv'
  usecase.run
end
