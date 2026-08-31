#!/usr/bin/env ruby
# encoding: utf-8
Version = '20260831-000000'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class SpliceWizApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'SpliceWiz'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Alternative_Splicing'
    @description =<<-EOS
    Detection and classification of differential alternative splicing
    (exon skipping, mutually-exclusive exons, alternative 5'/3' splice sites,
    intron retention, alternative first/last exons) from STAR BAM files, for a
    two-group comparison.<br/>
<a href='https://bioconductor.org/packages/release/bioc/html/SpliceWiz.html'>SpliceWiz</a><br/>
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
    @params['grouping'] = ''
    @params['sampleGroup'] = ''
    @params['sampleGroup', 'description'] = 'sampleGroup should be different from refGroup'
    @params['refGroup'] = ''
    @params['refGroup', 'description'] = 'refGroup should be different from sampleGroup'
    @params['aseMethod'] = ['limma', 'DESeq', 'edgeR']
    @params['aseMethod', 'description'] = 'differential ASE engine (ASE_limma / ASE_DESeq / ASE_edgeR)'
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['strandMode', 'description'] = "library strandedness; 'both' collates strand-agnostically, otherwise SpliceWiz auto-detects the direction"
    @params['IRmode'] = ['all', 'annotated', 'annotated_binary']
    @params['IRmode', 'description'] = 'intron-retention handling for the ASE test'
    @params['batch'] = ''
    @params['batch', 'description'] = 'optional Factor column to include as a batch covariate (~ batch + condition); leave empty for none'
    @params['FDR'] = '0.05'
    @params['FDR', 'description'] = 'adjusted-p cutoff for calling a significant event'
    @params['deltaPSI'] = '0.1'
    @params['deltaPSI', 'description'] = 'minimum |delta PSI| for calling a significant event'
    @params['topN'] = '20'
    @params['topN', 'description'] = 'number of top events to draw coverage plots for'
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
    command << run_RApp("EzAppSpliceWiz")
  end
end

if __FILE__ == $0

end
