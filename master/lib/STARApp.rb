#!/usr/bin/env ruby
# encoding: utf-8
Version = '20170420-093005'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class STARApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'STAR'
    @analysis_category = 'Map'
    @description =<<-EOS
    Ultafast spliced alignment<br/>
<a href='https://code.google.com/p/rna-star/'>manual</a><br/>
Noteworthy options:<ul>
<li> --outFilterMatchNmin 30 --outFilterMismatchNmax 5 --outFilterMismatchNoverLmax 0.05 --outFilterMultimapNmax 50</li>
<li> large numbers of contigs ask form more RAM. In that case the index must be built with a smaller --genomeChrBinNbits 18;
e.g. the Wheat assembly with 1Mio contigs and --genomeChrBinNbits 10 still takes 120GB of RAM during alignmnet</li>
</ul>
EOS
    @required_columns = ['Name','Read1','Species']
    @required_params = ['refBuild','paired', 'strandMode']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '30'
    @params['scratch'] = '100'
    @params['refBuild'] = ref_selector
    @params['paired'] = false
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['refFeatureFile'] = 'genes.gtf'
    @params['cmdOptions'] = '--outFilterType BySJout --outFilterMatchNmin 30 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.05 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignIntronMax 1000000 --alignMatesGapMax 1000000  --outFilterMultimapNmax 50 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --outSAMstrandField intronMotif'
    @params['getChimericJunctions'] = false
    @params['trimAdapter'] = true
    @params['trimLeft'] = 0
    @params['trimRight'] = 0
    @params['minTailQuality'] = 0
    @params['minAvgQuality'] = 10
    @params['minReadLength'] = 20
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Aligner/STAR", "Tools/samtools"]
  end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  end
  def set_default_parameters
    if dataset_has_column?('refBuild')
      @params['refBuild'] = @dataset[0]['refBuild']
    end
    @params['paired'] = dataset_has_column?('Read2')
    if dataset_has_column?('strandMode')
      @params['strandMode'] = @dataset[0]['strandMode']
    end
  end
  def next_dataset
     dataset = {
        'Name'=>@dataset['Name'], 
        'BAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam"), 
        'BAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam.bai"),
        'IGV Starter [Link]'=>File.join(@result_dir, "#{@dataset['Name']}-igv.jnlp"),
        'Species'=>@dataset['Species'],
        'refBuild'=>@params['refBuild'],
        'paired'=>@params['paired'],
        'refFeatureFile'=>@params['refFeatureFile'],
        'strandMode'=>@params['strandMode'],
        'Read Count'=>@dataset['Read Count'],
        'IGV Starter [File]'=>File.join(@result_dir, "#{@dataset['Name']}-igv.jnlp"),
        'IGV Session [File]'=>File.join(@result_dir, "#{@dataset['Name']}-igv.xml"),
        'PreprocessingLog [File]'=>File.join(@result_dir, "#{@dataset['Name']}_preprocessing.log"),
        'STARLog [File]'=>File.join(@result_dir, "#{@dataset['Name']}_STAR.log")
     }.merge(extract_column("Factor")).merge(extract_column("B-Fabric"))

     if @params['getChimericJunctions']
       dataset['Chimerics [File]'] = File.join(@result_dir, "#{@dataset['Name']}.chimeric")
     end
     dataset
  end
  def commands
    run_RApp("EzAppSTAR")
  end
end

if __FILE__ == $0
  run STARApp

end

