#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-095250'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class RSEMApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'RSEM'
    @analysis_category = 'Count'
    @description =<<-EOS
    Use bowtie alignments to transcript database and a posterior model to estimate isoform/gene abundances<br/>
<a href='http://deweylab.biostat.wisc.edu/rsem/rsem-calculate-expression.html'>manual/</a><br/>
Noteworthy is the option --bowtie-e which can be used to limit the sum of mismatching qualities for the alignments
EOS
    @required_columns = ['Name','Read1','Species']
    @required_params = ['paired', 'strandMode']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '100'
    @params['refBuild'] = ref_selector
    @params['paired'] = false
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['refFeatureFile'] = 'genes.gtf'
    @params['bowtie-e'] = '200'
    @params['bowtie-e', 'description'] = 'maximum sum of base qualities at mismatching positions'
    @params['cmdOptions'] = ' --calc-ci --sort-bam-by-read-name'
    @params['keepBam'] = false
    @params['keepBam', 'description'] = 'converts the transcript alignments into genome coordinates and reports them as a BAM file'
    @params['trimAdapter'] = true
    @params['trimLeft'] = 0
    @params['trimRight'] = 0
    @params['minTailQuality'] = 0
    @params['minAvgQuality'] = 20
    @params['specialOptions'] = ''
    @params['transcriptFasta'] = ''
    @params['transcriptFasta', 'description'] = 'give full path of transcript fasta file; in that case the build is ignored; if it comes from trinity assembly the gene-isoform associations will be extracted and used'
    @params['transcriptTypes'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA', 'long_noncoding', 'short_noncoding', 'pseudogene']
    @params['transcriptTypes', 'multi_selection'] = true
    @params['transcriptTypes', 'selected'] = 'protein_coding'
    @params['mail'] = ""
	  # Bowtie >=1.2.0 may return interleaving mates which trips RSEM (as of v1.3.0) as it expects
	  # each read to be followed by a mate.
	  # Also, as of v1.3.0, it only supports samtools v1.3.1
	  @modules = ["Tools/samtools/1.3.1", "Aligner/Bowtie/1.1.2", "Aligner/Bowtie2", "Aligner/STAR", "Aligner/RSEM", "QC/Flexbar", "QC/Trimmomatic", "Dev/R", "Tools/sambamba"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  end
  def set_default_parameters
    @params['paired'] = dataset_has_column?('Read2')
    if dataset_has_column?('strandMode')
      @params['strandMode'] = @dataset[0]['strandMode']
    end
  end
  def next_dataset
    if @params['keepBam']
      if @params['transcriptFasta'] == ''
        {'Name'=>@dataset['Name'],
         'Count [File]'=>File.join(@result_dir, "#{@dataset['Name']}.txt"),
         'BAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam"),
         'BAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam.bai"),
         'Species'=>@dataset['Species'],
         'refBuild'=>@params['refBuild'],
         'featureLevel'=>'isoform',
         'refFeatureFile'=>@params['refFeatureFile'],
         'strandMode'=>@params['strandMode'],
         'paired'=>@params['paired'],
         'Read Count'=>@dataset['Read Count'],
         'transcriptTypes'=>@params['transcriptTypes']
        }.merge(extract_columns(@inherit_tags))
      else
        {'Name'=>@dataset['Name'],
         'Count [File]'=>File.join(@result_dir, "#{@dataset['Name']}.txt"),
         'BAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam"),
         'BAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam.bai"),
         'Species'=>@dataset['Species'],
         'refBuild'=>@params['refBuild'],
         'featureLevel'=>'isoform',
         'refFeatureFile'=>@params['refFeatureFile'],
         'strandMode'=>@params['strandMode'],
         'paired'=>@params['paired'],
         'Read Count'=>@dataset['Read Count'],
         'transcriptTypes'=>'NA'
        }.merge(extract_columns(@inherit_tags))
      end
    else
      if @params['transcriptFasta'] == ''
        {'Name'=>@dataset['Name'],
         'Count [File]'=>File.join(@result_dir, "#{@dataset['Name']}.txt"),
         'Species'=>@dataset['Species'],
         'refBuild'=>@params['refBuild'],
         'featureLevel'=>'isoform',
         'refFeatureFile'=>@params['refFeatureFile'],
         'strandMode'=>@params['strandMode'],
         'paired'=>@params['paired'],
         'Read Count'=>@dataset['Read Count'],
         'transcriptTypes'=>@params['transcriptTypes'],
         'PreprocessingLog [File]'=>File.join(@result_dir, "#{@dataset['Name']}_preprocessing.log")
        }.merge(extract_columns(@inherit_tags))
      else
        {'Name'=>@dataset['Name'],
         'Count [File]'=>File.join(@result_dir, "#{@dataset['Name']}.txt"),
         'Species'=>@dataset['Species'],
         'refBuild'=>@params['refBuild'],
         'featureLevel'=>'isoform',
         'refFeatureFile'=>@params['refFeatureFile'],
         'strandMode'=>@params['strandMode'],
         'paired'=>@params['paired'],
         'Read Count'=>@dataset['Read Count'],
         'transcriptTypes'=>'NA',
         'PreprocessingLog [File]'=>File.join(@result_dir, "#{@dataset['Name']}_preprocessing.log")
        }.merge(extract_columns(@inherit_tags))
      end
    end
  end
  def commands
    run_RApp("EzAppRSEM")
  end
end

if __FILE__ == $0
  run RSEMApp

end
