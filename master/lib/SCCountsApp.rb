#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-095351'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class SCCountsApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'SCCountsApp'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
    Maps unmapped bam file with STAR and generates stats and expression counts with featureCounts
EOS
    @required_columns = ['Name','Read1','Species', 'CellDataset']
    @required_params = ['refBuild','paired']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '50'
    @params['scratch'] = '500'
    @params['refBuild'] = ref_selector
    @params['paired'] = false
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['refFeatureFile'] = 'genes.gtf'
    @params['mapMethod'] = ['STAR', 'bowtie', 'bowtie2', 'tophat', 'bwa-mem']
    @params['mapOptions'] = '--outFilterType BySJout --outFilterMatchNmin 30 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.05 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignIntronMax 1000000 --alignMatesGapMax 1000000  --outFilterMultimapNmax 50 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --outSAMstrandField intronMotif'
    @params['getChimericJunctions'] = false
    @params['trimAdapter'] = true
    @params['trimLeft'] = 0
    @params['trimRight'] = 0
    @params['minTailQuality'] = 10
    @params['minAvgQuality'] = 10
    @params['minReadLength'] = 20
    @params['featureLevel'] = 'gene'
    @params['gtfFeatureType'] = 'exon'
    @params['allowMultiOverlap'] = true
    @params['allowMultiOverlap', 'description'] = "count alignments that fall in a region where multipe features are annotated"
    @params['countPrimaryAlignmentsOnly'] = true
    @params['minFeatureOverlap'] = 10
    @params['minFeatureOverlap', 'description'] = "minimum overlap of a read with a transcript feature"
    @params['minMapQuality'] = 1
    @params['keepMultiHits'] = true
    @params['transcriptTypes'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA', 'long_noncoding', 'short_noncoding', 'pseudogene']
    @params['transcriptTypes', 'multi_selection'] = true
    @params['transcriptTypes', 'selected'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA']
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Dev/jdk", "Aligner/STAR", "Tools/samtools", "Aligner/BWA", "Aligner/Bowtie", "Aligner/Bowtie2", "Aligner/TopHat", "QC/Trimmomatic", "QC/Flexbar", "Tools/Picard", "Dev/Python2", "Dev/R", "Tools/sambamba"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'paired'=>@params['paired'],
     'featureLevel'=>@params['featureLevel'],
     'refFeatureFile'=>@params['refFeatureFile'],
     'transcriptTypes'=>@params['transcriptTypes'],
     'CellDataset [File]'=>File.join(@result_dir, "#{@dataset['Name']}-dataset.tsv"),
     'CountMatrix [File]'=>File.join(@result_dir, "#{@dataset['Name']}-counts.txt"),
     'Stats [File]'=>File.join(@result_dir, "#{@dataset['Name']}-stats.txt"),
     'CellCyclePhase [File]'=>File.join(@result_dir, "#{@dataset['Name']}-CellCyclePhase.txt"),
     'BAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam"),
     'BAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam.bai"),
     'PreprocessingLog [File]'=>File.join(@result_dir, "#{@dataset['Name']}_preprocessing.log"),
     'STARLog [File]'=>File.join(@result_dir, "#{@dataset['Name']}_STAR.log")
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppSCCounts")
  end
end
