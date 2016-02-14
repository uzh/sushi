#!/usr/bin/env ruby
# encoding: utf-8
Version = '20160215-004949'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class MACS2App < SushiFabric::SushiApp
  def initialize
    super
    @name = 'MACS2'
    @analysis_category = 'Peaks'
    @description =<<-EOS
Capturing the influence of genome complexity to evaluate the significance of enriched ChIP regions<br/>
<a href='http://liulab.dfci.harvard.edu/MACS/00README.html'>http://liulab.dfci.harvard.edu/MACS/00README.html</a>
    EOS
    @required_columns = ['Name','BAM','BAI', 'refBuild','Control']
    @required_params = ['refBuild','paired']
    # optional params
    @params['cores'] = '1'
    @params['ram'] = '16'
    @params['scratch'] = '100'
    @params['refBuild'] = ref_selector
    @params['paired'] = false
    @params['useControl'] = true
    @params['refFeatureFile'] = 'genes.gtf'
    @params['cmdOptions'] = '--nomodel --extsize 147 -g hs --bw 200'
    @params['specialOptions'] = ''
    @params['mail'] = ''
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('refFeatureFile')
      @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    end
    if dataset_has_column?('paired')
      @params['paired'] = @dataset[0]['paired']
    end                               
 end

  def next_dataset
    bw_link = File.join(@result_dir, "#{@dataset['Name']}.bw")
    peakfile_link = File.join(@result_dir, "#{@dataset['Name']}_peaks.xls")
    bedfile_link = File.join(@result_dir, "#{@dataset['Name']}_peaks.bed")
    peakseq_link = File.join(@result_dir, "#{@dataset['Name']}_peaks.fa")

    {'Name'=>@dataset['Name'], 
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'refFeatureFile'=>@params['refFeatureFile'],
     'paired'=>@params['paired'],
     'CalledPeaks [File]'=>peakfile_link,
     'BED [File]'=>bedfile_link,
     'PeakSequences [File]'=>peakseq_link,
     'BigWigFile [File]'=>bw_link
    }.merge(extract_column("Factor")).merge(extract_column("B-Fabric")).merge(extract_column("ChIP"))
  end
  def commands
    run_RApp("EzAppMacs2")
  end
end

if __FILE__ == $0

end

