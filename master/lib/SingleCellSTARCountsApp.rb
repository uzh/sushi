#!/usr/bin/env ruby
# encoding: utf-8
Version = '20150507-090225'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class SingleCellSTARCountsApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'SingleCellSTARCountsApp'
    @analysis_category = 'Count'
    @description =<<-EOS
    Maps all read files found in a folder with STAR and generates stats and expression counts with featureCounts
EOS
    @required_columns = ['Name','ReadFolder','Species']
    @required_params = ['refBuild','paired']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '300'
    @params['refBuild'] = ref_selector
    @params['paired'] = false
#    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['refFeatureFile'] = 'genes.gtf'
    @params['cmdOptions'] = '--genomeLoad LoadAndKeep --outFilterType BySJout --outFilterMatchNmin 30 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.05 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignIntronMax 1000000 --alignMatesGapMax 1000000  --outFilterMultimapNmax 50 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --outSAMstrandField intronMotif'
    @params['getChimericJunctions'] = false
    @params['trimAdapter'] = true
    @params['trimLeft'] = 3
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
    @params['minMapQuality'] = 10
    @params['keepMultiHits'] = true
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Aligner/STAR", "Tools/samtools"]

  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'CountFolder [File]'=>File.join(@result_dir, "#{@dataset['Name']}"),
     'CountDataset [Link]'=>File.join(@result_dir, "#{@dataset['Name']}/#{@dataset['Name']}-dataset.tsv"),
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'paired'=>@params['paired'],
     'refFeatureFile'=>@params['refFeatureFile'],
      }.merge(extract_column("Factor")).merge(extract_column("B-Fabric"))
  end
  def commands
    command = "#{R_COMMAND} --vanilla --slave<<  EOT\n"
    command << "EZ_GLOBAL_VARIABLES <<- '#{EZ_GLOBAL_VARIABLES}'\n"
    command << "library(ezRun)\n"
    command << "param = list()\n"
    param = @params
    param.keys.each do |key|
      command << "param[['#{key}']] = '#{param[key]}'\n"
    end
    command << "param[['dataRoot']] = '#{@gstore_dir}'\n"
    command << "param[['resultDir']] = '#{@result_dir}'\n"
    command << "param[['writeIgvSessionLink']] = FALSE\n"
    command << "output = list()\n"
    output = next_dataset
    output.keys.each do |key|
      command << "output[['#{key}']] = '#{output[key]}'\n"
    end
    command << "input = list()\n" 
    input = @dataset
    input.keys.each do |key|
      command << "input[['#{key}']] = '#{input[key]}'\n" 
    end
    command << "countDir = file.path(getwd(), basename(output[['CountFolder [File]']]))\n"
    command << "dir.create(countDir)\n"
    command << "param = ezParam(param)\n"
    command << "dsDir = input[['ReadFolder']]\n"
    command << "old = getwd()\n"
    command << "setwd(param[['dataRoot']])\n"
    command << "if (param[['paired']]){\n"
    command << "  meta = makeMinimalPairedEndReadDataset(fqDir = dsDir, species='Mus_musculus', strandMode = 'both', readCount = NULL, readTypeSuffix= c('_1.fastq.gz', '_2.fastq.gz'))\n"
    command << "} else {\n"
    command << "  meta = makeMinimalSingleEndReadDataset(fqDir = dsDir, species='Mus_musculus', strandMode = 'both', readCount = NULL)\n"
    command << "}\n"
    command << "rownames(meta) = basename(rownames(meta))\n"
    command << "setwd(old)\n"
    command << "readDs = EzDataset(meta=meta, dataRoot=param[['dataRoot']])\n"
    command << "meta[['BAM [File]']] = paste0(rownames(meta), '.bam')\n"
    command << "meta[['BAI [File]']] = paste0(rownames(meta), '.bam.bai')\n"
    command << "bamDs = EzDataset(meta=meta)\n"
    command << "meta[['Count [File]']] = paste0(rownames(meta), '.txt')\n"
    command << "meta[['Stats [File]']] = paste0(rownames(meta), '-stats.txt')\n"
    command << "ezWrite.table(meta, file=file.path(countDir, basename(output[['CountDataset [Link]']])), head='Name')\n"
    command << "countDs = EzDataset(meta=meta)\n"
    command << "mailAddress = param[['mail']]\n"
    command << "param[['mail']]  = ''\n"
    command << "# copy the index locally\n"
    command << "#refDir = sub('.gtf$', '_STARIndex', param$ezRef['refFeatureFile'])\n"
    command << "#ezSystem(paste('cp -r', refDir, basename(refDir)))\n"
    command << "#param$ezRef['refIndex'] = file.path(getwd(), basename(refDir))\n"
    command << "for (nm in readDs\\$getNames()){\n"
    command << "  message(nm)\n" 
    command << "  setwdNew(nm)\n"
    command << "  readI = readDs\\$subset(nm)\n"
    command << "  bamI = bamDs\\$subset(nm)\n"
    command << "  countI = countDs\\$subset(nm)\n"
    command << "  EzAppSTAR\\$new()\\$run(input=readI, output=bamI, param=param)\n"
    command << "  starLogfile = paste0(countDir, '/', nm, '-STAR.log')\n"
    command << "  file.rename('Log.final.out', to = starLogfile)\n"
    command << "  for (strandMode in c('sense', 'antisense', 'both')){\n"
    command << "    param\\$strandMode = strandMode\n"
    command << "    countI\\$setColumn('Count', paste0(nm, '-', param\\$strandMode, '.txt'))\n"
    command << "    countI\\$setColumn('Stats', paste0(nm, '-', param\\$strandMode, '-stats.txt'))\n"
    command << "    EzAppFeatureCounts\\$new()\\$run(input=bamI, output=countI, param=param)\n"
    command << "  }\n"
    command << "  file.copy(list.files('.', paste0(nm, '.*txt')), countDir)\n"
    command << "  setwd('..')\n"
    command << "  unlink(nm, recursive=TRUE, force=TRUE)\n"
    command << "}\n"
    command << "refDir = sub('.gtf$', '_STARIndex', param$ezRef['refFeatureFile'])\n"
    command << "ezSystem(paste('STAR', '--genomeDir', refDir, '--genomeLoad Remove'))\n"
    command << "if (ezValidMail(mailAddress)){\n"
    command << "  ezMail(subject=paste(basename(dsDir), ' -- single cell counts done'), text=paste(basename(dsDir), ' -- single cell counts done'), to=mailAddress)\n"
    command << "}\n"
    command <<  "EOT\n"
    command
  end
end

