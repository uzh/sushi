#!/usr/bin/env ruby
# encoding: utf-8
Version = '20170926-161518'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class DiscovarDenovoApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'DiscovarDenovo'
    @analysis_category = 'Assemble'
    @params['process_mode'] = 'DATASET'
    @description =<<-EOS
DISCOVAR de novo is a large (and small) de novo genome assembler. It quickly generates highly accurate and complete assemblies using the same single library data as used by DISCOVAR. It currently doesn’t support variant calling – for that, please use DISCOVAR instead.
Refer to <a href='https://software.broadinstitute.org/software/discovar/blog/'>https://software.broadinstitute.org/software/discovar/blog/</a>
    EOS
    @required_columns = ['Name','Read1']
    @required_params = []
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '16'
    @params['scratch'] = '100'
    @modules = ["Assembly/DiscovarDenovo"]
    #@inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def next_dataset
    {'Name'=>"Assembly", 
     'Results [File]'=>File.join(@result_dir, "discovar_out"),	
    }
  end
  def commands
    command = ""
    command << "mkdir discovar_out\n"
    file_cols = get_columns_with_tag("File")
    reads = file_cols.map{|data_set| data_set.values.map{|file| File.join(@gstore_dir, file)}}.flatten.join(",")
    command << "DiscovarDeNovo READS=#{reads} OUT_DIR=discovar_out NUM_THREADS=#{@params['cores']} MAX_MEM_GB=#{@params['ram']}\n"
    command
  end
end

if __FILE__ == $0

end

