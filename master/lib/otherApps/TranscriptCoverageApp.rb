#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-095715'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class TranscriptCoverageApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'TranscriptCoverage'
    @analysis_category = 'Profiles'
 @description =<<-EOS
    Compute the base-by-base coverage from transript alignment files.
EOS
    @required_columns = ['Name','trBAM','trBAI', 'refBuild']
    @required_params = ['refBuild','paired', 'strandMode']
    # optional params
    @params['cores'] = '2'
    @params['ram'] = '30'
    @params['scratch'] = '10'
    @params['refBuild'] = ref_selector
    @params['paired'] = false
    @params['strandMode'] = ['both', 'sense', 'antisense']
    @params['refFeatureFile'] = 'genes.gtf'
    @params['minReadLength'] = 28
    @params['maxReadLength'] = 32
    @params['getCoverageByReadLength'] = true
    @params['coverageType'] = ['readStart', 'fullRead']
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
    @modules = ["Dev/R"]
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
     'Profiles [File]'=>File.join(@result_dir, "#{@dataset['Name']}-profiles.RData"),
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'refFeatureFile'=>@params['refFeatureFile'],
     'featureLevel'=>'isoform',
     'strandMode'=>@params['strandMode'],
     'paired'=>@params['paired'],
     'Read Count'=>@dataset['Read Count'],
     'coverageType'=>@params['coverageType']
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppTranscriptCoverage")
  end
end

if __FILE__ == $0

end

