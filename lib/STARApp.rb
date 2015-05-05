#!/usr/bin/env ruby
# encoding: utf-8
Version = '20150226-111433'

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
    @required_params = ['build','paired', 'strandMode']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '40'
    @params['scratch'] = '100'
    @params['build'] = ref_selector
    @params['paired'] = false
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['featureFile'] = 'genes.gtf'
    @params['cmdOptions'] = '--outFilterType BySJout --outFilterMatchNmin 30 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.05 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignIntronMax 1000000 --alignMatesGapMax 1000000  --outFilterMultimapNmax 50 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --outSAMstrandField intronMotif'
    @params['getChimericReads'] = false
    @params['trimAdapter'] = false
    @params['trimLeft'] = 0
    @params['trimRight'] = 0
    @params['minTailQuality'] = 0
    @params['specialOptions'] = ''
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
     if @params['getChimericReads']
       {'Name'=>@dataset['Name'], 
        'BAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam"), 
        'BAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam.bai"),
        'IGV Starter [Link]'=>File.join(@result_dir, "#{@dataset['Name']}-igv.jnlp"),
        'Species'=>@dataset['Species'],
        'build'=>@params['build'],
        'paired'=>@params['paired'],
        'featureFile'=>@params['featureFile'],
        'strandMode'=>@params['strandMode'],
        'Read Count'=>@dataset['Read Count'],
        'Chimerics [File]'=>File.join(@result_dir, "#{@dataset['Name']}.chimeric"),
        'IGV Starter [File]'=>File.join(@result_dir, "#{@dataset['Name']}-igv.jnlp"),
        'IGV Session [File]'=>File.join(@result_dir, "#{@dataset['Name']}-igv.xml")
      }.merge(extract_column("Factor")).merge(extract_column("B-Fabric"))
     else 
       {'Name'=>@dataset['Name'],
        'BAM [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam"),
        'BAI [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bam.bai"),
        'IGV Starter [Link]'=>File.join(@result_dir, "#{@dataset['Name']}-igv.jnlp"),
        'Species'=>@dataset['Species'],
        'build'=>@params['build'],
        'paired'=>@params['paired'],
        'featureFile'=>@params['featureFile'],
        'strandMode'=>@params['strandMode'],
        'Read Count'=>@dataset['Read Count'],
        'IGV Starter [File]'=>File.join(@result_dir, "#{@dataset['Name']}-igv.jnlp"),
        'IGV Session [File]'=>File.join(@result_dir, "#{@dataset['Name']}-igv.xml")
      }.merge(extract_column("Factor")).merge(extract_column("B-Fabric"))
     end
  end
  def commands
    run_RApp("mapSTARApp")
  end
end

if __FILE__ == $0
  run STARApp

end

