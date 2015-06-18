#!/usr/bin/env ruby
# encoding: utf-8
Version = '20150226-110907'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class FeatureCountsApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'FeatureCounts'
    @analysis_category = 'Count'
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
    @params['allowMultiOverlap' = true
    @params['allowMultiOverlap', 'description'] = "count alignments that fall in a region where multipe features are annotated"
    @params['countPrimaryAlignmentsOnly'] = true
    @params['minFeatureOverlap'] = 10
    @params['minFeatureOverlap', 'description'] = "minimum overlap of a read with a transcript feature"
    @params['countTrimmedTranscripts'] = false
    @params['trimmedMaxLength'] = 200
    @params['minMapQuality'] = 10
    @params['keepMultiHits'] = true
    @params['specialOptions'] = ''
    @params['mail'] = ""
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
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'featureLevel'=>@params['featureLevel'],
     'refFeatureFile'=>@params['refFeatureFile'],
     'strandMode'=>@params['strandMode'],
     'paired'=>@params['paired'],
     'Read Count'=>@dataset['Read Count']
    }.merge(extract_column("Factor")).merge(extract_column("B-Fabric"))
  end
  def commands
    run_RApp("ezAppFeatureCounts")
  end
end

if __FILE__ == $0

end

