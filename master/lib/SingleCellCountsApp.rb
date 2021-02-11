#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-095406'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class SingleCellCountsApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'SingleCellCountsApp'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
    Maps all read files specified by the dataset file generates stats and expression counts with featureCounts
EOS
    @required_columns = ['Name','ReadDataset','Species']
    @required_params = ['refBuild','paired']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '300'
    @params['refBuild'] = ref_selector
    @params['paired'] = false
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['refFeatureFile'] = 'genes.gtf'
    @params['mapMethod'] = ['STAR', 'bowtie', 'bowtie2', 'tophat', 'bwa-mem']
    @params['mapOptions'] = '--genomeLoad LoadAndKeep --outFilterType BySJout --outFilterMatchNmin 30 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.05 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignIntronMax 1000000 --alignMatesGapMax 1000000  --outFilterMultimapNmax 50 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --outSAMstrandField intronMotif'
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
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Dev/jdk", "Aligner/STAR", "Tools/samtools", "Tools/Picard", "Dev/Python2", "Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'CountDataset [File]'=>File.join(@result_dir, "#{@dataset['Name']}-dataset.tsv"),
     'CountMatrix [File]'=>File.join(@result_dir, "#{@dataset['Name']}-counts.txt"),
     'CountFolder [File]'=>File.join(@result_dir, "#{@dataset['Name']}-Counts"),
     'CellCyclePhase [File]'=>File.join(@result_dir, "#{@dataset['Name']}-CellCyclePhase.txt"),
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'paired'=>@params['paired'],
     'refFeatureFile'=>@params['refFeatureFile'],
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppSingleCellCounts")
  end
end
