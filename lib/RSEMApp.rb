#!/usr/bin/env ruby
# encoding: utf-8
Version = '20150226-111524'

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
    @params['cmdOptions'] = ' --calc-ci '
    @params['keepBam'] = false
    @params['trimAdapter'] = true
    @params['trimLeft'] = 0
    @params['trimRight'] = 0
    @params['minTailQuality'] = 0
    @params['minAvgQuality'] = 20
    @params['specialOptions'] = ''
    @params['trinityFasta'] = ''
    @params['trinityFasta', 'description'] = 'give full path of trinity assembled fasta output file; in that case the build is ignored'
    @params['mail'] = ""
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
       'Read Count'=>@dataset['Read Count']
      }.merge(extract_column("Factor")).merge(extract_column("B-Fabric"))
    else 
      {'Name'=>@dataset['Name'],
       'Count [File]'=>File.join(@result_dir, "#{@dataset['Name']}.txt"),
       'Species'=>@dataset['Species'],
       'refBuild'=>@params['refBuild'],
       'featureLevel'=>'isoform',
       'refFeatureFile'=>@params['refFeatureFile'],
       'strandMode'=>@params['strandMode'],
       'paired'=>@params['paired'],
       'Read Count'=>@dataset['Read Count']
      }.merge(extract_column("Factor")).merge(extract_column("B-Fabric"))
    end
  end
  def commands
    run_RApp("EzAppRSEM")
  end
end

if __FILE__ == $0
  run RSEMApp

end

