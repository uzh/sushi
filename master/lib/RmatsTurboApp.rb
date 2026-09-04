#!/usr/bin/env ruby
# encoding: utf-8
Version = '20260902-000000'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class RmatsTurboApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'rMATS-turbo'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Alternative_Splicing'
    @description =<<-EOS
    Detection and quantification of differential alternative splicing
    (exon skipping, mutually-exclusive exons, alternative 5'/3' splice sites,
    intron retention) from STAR BAM files, for a two-group comparison, with
    rMATS-turbo (replicate MATS).<br/>
<a href='https://github.com/Xinglab/rmats-turbo'>rMATS-turbo</a><br/>
    EOS
    @required_columns = ['Name','BAM','BAI','Species','refBuild','refFeatureFile']
    @required_params = ['grouping','sampleGroup','refGroup','refBuild']
    # optional params
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '40'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '100'
    @params['scratch', "context"] = "slurm"
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "reference genome assembly"
    @params['refFeatureFile'] = 'genes.gtf'
    @params['transcriptTypes'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA', 'long_noncoding', 'short_noncoding', 'pseudogene']
    @params['transcriptTypes', 'multi_selection'] = true
    @params['transcriptTypes', 'selected'] = 'protein_coding'
    @params['transcriptTypes', 'description'] = 'restrict the reference GTF to these transcript biotypes before event discovery (fewer non-coding events)'
    @params['grouping'] = ''
    @params['sampleGroup'] = ''
    @params['sampleGroup', 'description'] = 'sampleGroup should be different from refGroup'
    @params['refGroup'] = ''
    @params['refGroup', 'description'] = 'refGroup should be different from sampleGroup'
    @params['readType'] = ['auto', 'paired', 'single']
    @params['readType', 'description'] = "rMATS -t; 'auto' reads the dataset 'paired' column"
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['strandMode', 'description'] = 'library strandedness -> rMATS --libType (both=fr-unstranded, sense=fr-secondstrand, antisense=fr-firststrand)'
    @params['readLength'] = ''
    @params['readLength', 'description'] = 'RNA-seq read length; leave blank to auto-detect from the BAMs (run with --variable-read-length)'
    @params['statModel'] = ['default', 'paired', 'darts']
    @params['statModel', 'description'] = 'rMATS statistical model: default (unpaired) / paired (--paired-stats) / darts (--darts-model)'
    @params['novelSS'] = ['false', 'true']
    @params['novelSS', 'description'] = 'detect novel (unannotated) splice sites (--novelSS)'
    @params['cstat'] = '0.0001'
    @params['cstat', 'description'] = 'rMATS --cstat: cutoff splicing difference for the null hypothesis test (0 <= c < 1)'
    @params['minCoverage'] = '10'
    @params['minCoverage', 'description'] = 'coverage filter: minimum average junction reads (inclusion+skipping, IJC+SJC) required in BOTH groups; rMATS is anti-conservative on low-coverage junctions, so events below this are greyed in the volcano and excluded from the FDR (multiple-testing) correction'
    @params['FDR'] = '0.05'
    @params['FDR', 'description'] = 'adjusted-p cutoff for calling a significant event (report-side)'
    @params['deltaPSI'] = '0.1'
    @params['deltaPSI', 'description'] = 'minimum |IncLevelDifference| for calling a significant event (report-side)'
    @params['topN'] = '20'
    @params['topN', 'description'] = 'number of top events to draw PSI boxplots for'
    @params['grouping2'] = ''
    @params['grouping2', 'description'] =  'optional secondary co-variate (reserved; rMATS models a single two-group factor)'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @params['Rversion'] = ["Dev/R/4.6.0"]
    @modules = ["Tools/samtools"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def preprocess
    @name = "#{@name}_#{@params['sampleGroup']}--over--#{@params['refGroup']}"
  end
  def next_dataset
    @comparison = "#{@params['sampleGroup']}--over--#{@params['refGroup']}"
    @params['comparison'] = @comparison
    @params['name'] = @comparison
    report_file = File.join(@result_dir, "#{@params['comparison']}")
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@comparison,
     'Species'=>(dataset = @dataset.first and dataset['Species']),
     'refBuild'=>@params['refBuild'],
     'Static Report [Link]'=>report_link,
     'Live Report [Link]'=>"http://fgcz-shiny.uzh.ch/exploreRmatsTurbo?data=#{report_file}",
     'Report [File]'=>report_file,
    }.merge(extract_columns(@inherit_tags))
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('refFeatureFile')
      @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    end
    if dataset_has_column?('strandMode')
      @params['strandMode'] = @dataset[0]['strandMode']
    end
  end
  def commands
    command = "module load #{@params["Rversion"]}\n"
    command << run_RApp("EzAppRmatsTurbo")
  end
end

if __FILE__ == $0

end
