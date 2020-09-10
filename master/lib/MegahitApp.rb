#!/usr/bin/env ruby
# encoding: utf-8
Version = '20181127-155223'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class MegahitApp < SushiFabric::SushiApp
def initialize
super
@name = 'Megahit'
@analysis_category = 'Metagenomics'
@description =<<-EOS
Denovo metagenomics assembly with Megahit. 
<a href='https://github.com/voutcn/megahit'>https://github.com/voutcn/megahit</a>
  EOS
@required_columns = ['Name', 'Read1']
@required_params = ['kmerMin','kmerMax','kmerStep','noMercy','kmin1pass','minCount']
@params['cores'] = '1'
@params['ram'] = '8'
@params['scratch'] = '10'
@params['kmerMin'] = '31'
@params['kmerMin', 'description'] = 'Minimum k-mer value for the assembly.'
@params['kmerMax'] = '101'
@params['kmerMax', 'description'] = 'Maximum k-mer value for the assembly.'
@params['kmerStep'] = '10'
@params['kmerStep', 'description'] = 'Step value to move from kmerMin to kmerMax.'
@params['noMercy'] = false
@params['noMercy', 'description'] = 'Recommended for assemblying metagenomes. Enable it if assemblying isolate with coverage above 30x.'
@params['kmin1pass'] = false
@params['kmin1pass', 'description'] = 'Recommended for generic metagenomes. Enable it for ultra-complex samples, e.g., soil.'
@params['minCount'] = '2'
@params['minCount', 'description'] = 'Recommended for assemblying metagenomes. Increase it to 3 or 4 if assemblying isolate with coverage above 40x.'
@params['mail'] = ""
@inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
@modules = ["Dev/R","Assembly/megahit"]
end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  __END__
  def set_default_parameters
     @params['paired'] = dataset_has_column?('Read2')
  end
def next_dataset
     {'Name'=>@dataset['Name'],
     'contigFile [File]' => File.join(@result_dir, "#{@dataset['Name']}.contigs.fasta"),
}
end
def commands
run_RApp("EzAppMegahit")
end
end

if __FILE__ == $0
end
