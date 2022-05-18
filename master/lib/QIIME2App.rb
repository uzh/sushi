#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class QIIME2App < SushiFabric::SushiApp
  def initialize
    super
    @name = 'QIIME2App'
    @analysis_category = 'Metagenomics'
    @description =<<-EOS
    Data processing with QIIME2. For short reads/Illumina data only.
    <a href='https://qiime2.org/'>QIIME2 main website with tutorials and manuals.</a>

EOS
    @params['process_mode'] = 'DATASET'
    @required_columns = ['Name', 'Read1']
    @required_params = ['paired','trim_left','truncate_len','sampling_depth', 'group']
    @params['cores'] = '2'
    @params['ram'] = '7'
    @params['scratch'] = '10'
    @params['paired'] = false
    @params['paired', 'description'] = 'whether the reads are paired end; if false then only Read1 is considered even if Read2 is available.'
    @params['trim_left'] = '0'
    @params['trim_left', 'description'] = 'Position at which sequences should be trimmed due to low quality. Default:0, assuming good quality reads.'
    @params['truncate_len'] = '150'
    @params['truncate_len', 'description'] = 'Position at which sequences should be truncated due to decrease in quality. Assuming good quality reads and 150bp illumina sequencing.'
    @params['sampling_depth'] = '1000'
    @params['sampling_depth', 'description'] = 'Total frequency that each sample should be rarefied to prior to computing diversity metrics.'
    @params['group'] = true
    @params['group', 'description'] = 'There needs to be a group assignment column. Ensure the column name is in the format "NAME [Factor]" and is placed as a column between Read1 and Read2'
    @params['name'] = 'QIIME2'
    @params['mail'] = ""
    @inherit_tags = ['B-Fabric']
    @modules = ['Dev/R']
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
     nds = {'Name'=>@params['name']}
     nds['Metadata [File]'] = File.join(@result_dir, 'sample_metadata.tsv')
     nds['Demux Report [Link]'] = File.join(@result_dir, 'demux_seqs.qzv.zip.folder/data/index.html')
     nds['Feature Table Html [Link]'] = File.join(@result_dir, 'table.qzv.zip.folder/data/index.html')
     nds['Rep Seqs Report [Link]'] = File.join(@result_dir, 'dada2_rep_set.qzv.zip.folder/data/index.html')
     nds['Denoising stats Html [Link]'] = File.join(@result_dir, 'dada2_denoising_stats.qzv.zip.folder/data/index.html')
     nds
  end
  def commands
     run_RApp("EzAppQIIME2", conda_env: "qiime2-2021.11")
  end
end

if __FILE__ == $0

end
