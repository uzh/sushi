#!/usr/bin/env ruby
# encoding: utf-8
Version = '20170907-111214'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class CountOverlapsApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'CountOverlaps'
    @analysis_category = 'Count'
    @required_columns = ['Name','BAM','BAI', 'refBuild']
    @required_params = ['refBuild','paired', 'strandMode']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '15'
    @params['scratch'] = '100'
    @params['refBuild'] = ref_selector
    @params['paired'] = false
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['refFeatureFile'] = 'genes.gtf'
    @params['featureLevel'] = 'gene'
    @params['countNonredundant'] = true
    @params['countNonredundant', 'description'] = "downweights alignments by the number of different genomic alignments"
    @params['minFeatureOverlap'] = 10
    @params['minFeatureOverlap', 'description'] = "minimum overlap of a read with a transcript feature"
    @params['countTrimmedTranscripts'] = false
    @params['trimmedMaxLength'] = 200
    @params['minMapQuality'] = 10
    @params['keepMultiHits'] = true
    @params['cmdOptions'] = ''
    @params['specialOptions'] = ''
    @params['mail'] = ""
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
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'featureLevel'=>@params['featureLevel'],
     'refFeatureFile'=>@params['refFeatureFile'],
     'strandMode'=>@params['strandMode'],
     'paired'=>@params['paired'],
     'Read Count'=>@dataset['Read Count']
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppCountOverlaps")
  end
end

if __FILE__ == $0

end

