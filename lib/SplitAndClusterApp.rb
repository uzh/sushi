#!/usr/bin/env ruby
# encoding: utf-8
Version = '20150710-112224'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class SplitAndClusterApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'SplitAndCluster'
    @analysis_category = 'Special'
    @description =<<-EOS
    Joins the paired end files, and then splits the joined file by internal forward and reverse primers. 
    Finally the resulting file is clustered with CD-Hit.<br>
    Use for p1871
    EOS
    @required_columns = ['Name','Read1', 'Read2']
    @required_params = ["forwardPrimerFile", "reversePrimerFile"]
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '16'
    @params['scratch'] = '100'
    @params['paired'] = true
    @params['paired', 'description'] = 'either the reads are paired-ends or single-end'
    @params['maxMismatch'] = '0'
    @params['maxMismatch', 'description'] = 'max number of mismatches in primer identification'
    @params['forwardPrimerFile'] = '/srv/gstore/projects/p1871/forwardPrimers_20160215.fasta'
    @params['reversePrimerFile'] = '/srv/gstore/projects/p1871/reversePrimers_20160215.fasta'
  end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  end
  def set_default_parameters
    @params['paired'] = dataset_has_column?('Read2')
  end
  def next_dataset
    nds = {'Name'=>@dataset['Name']}
    nds['Clustered [File]'] = File.join(@result_dir, "#{@dataset['Name']}_result")
    nds.merge(extract_column("Factor")).merge(extract_column("B-Fabric"))
    nds
  end
  def commands
    run_RApp("EzAppSplitAndCluster")
  end
end

if __FILE__ == $0

end

