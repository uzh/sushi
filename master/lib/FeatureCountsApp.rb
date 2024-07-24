#!/usr/bin/env ruby
# encoding: utf-8
Version = '20180209-163825'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class FeatureCountsApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'FeatureCounts'
    @analysis_category = 'Count'
 @description =<<-EOS
    Multi-purpose read counting with Rsubread::featureCounts<br/>
<a href='http://www.bioconductor.org/packages/release/bioc/manuals/Rsubread/man/Rsubread.pdf'>manual</a><br/>
Consider default values as suggestions. They are subject to change!
EOS
    @required_columns = ['Name','BAM','BAI', 'refBuild']
    @required_params = ['refBuild','paired', 'strandMode']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '20'
    @params['scratch'] = '10'
    @params['refBuild'] = ref_selector
    @params['paired'] = false
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['refFeatureFile'] = 'genes.gtf'
    @params['featureLevel'] = 'gene'
    @params['gtfFeatureType'] = 'exon'
    @params['gtfFeatureType', 'description'] = "which atomic features of the gtf should be used to define the meta-features; see featureLevel"
    @params['allowMultiOverlap'] = true
    @params['allowMultiOverlap', 'description'] = "count alignments that fall in a region where multipe features are annotated"
    @params['countPrimaryAlignmentsOnly'] = true
    @params['minFeatureOverlap'] = 10
    @params['minFeatureOverlap', 'description'] = "minimum overlap of a read with a transcript feature"
    @params['minMapQuality'] = 10
    @params['keepMultiHits'] = true
    @params['ignoreDup'] = false
    @params['transcriptTypes'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA', 'long_noncoding', 'short_noncoding', 'pseudogene']
    @params['transcriptTypes', 'multi_selection'] = true
    @params['transcriptTypes', 'selected'] = 'protein_coding'
    @params['secondRef'] = ''
    @params['secondRef', 'description'] = "extra DNA/RNA sequences to use for alignment; needs to point to a file on FGCZ servers; if the .fasta file has a corresponding .gtf file, this file needs to have the same base name, e.g. a file 'foo.fa' in folder /path/to/file/ requires a file 'foo.gtf' in the same folder in order for the gtf file to be used; ask for upload sushi@fgcz.ethz.ch."
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Tools/samtools", "Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('refFeatureFile')
      @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    end
    if dataset_has_column?('paired')
      @params['paired'] = @dataset[0]['paired']
    end
    if dataset_has_column?('strandMode')
      @params['strandMode'] = @dataset[0]['strandMode']
    end
  end

  def next_dataset
    {'Name'=>@dataset['Name'],
     'Count [File]'=>File.join(@result_dir, "#{@dataset['Name']}.txt"),
     'Stats [File]'=>File.join(@result_dir, "#{@dataset['Name']}-stats.txt"),
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'featureLevel'=>@params['featureLevel'],
     'refFeatureFile'=>@params['refFeatureFile'],
     'strandMode'=>@params['strandMode'],
     'paired'=>@params['paired'],
     'Read Count'=>@dataset['Read Count'],
     'transcriptTypes'=>@params['transcriptTypes']
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppFeatureCounts")
  end
end

if __FILE__ == $0

end
