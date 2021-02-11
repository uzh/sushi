#!/usr/bin/env ruby
# encoding: utf-8
Version = '20160425-070243'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class FastqcExample2 < SushiFabric::SushiApp
  def initialize
    super
    @name = 'FastqcExample2'
    @analysis_category = 'Map'
    @description =<<-EOS

    EOS
    @required_columns = ['Name','Read1']
    @required_params = ['paired']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '15'
    @params['scratch'] = '100'
    @params['library_type'] = ['fr-unstranded', 'fr-firststrand', 'fr-secondstrand']
    @params['paired'] = false
    @params['paired', 'description'] = 'whether the reads are paired-end or single-end'
  end
  def set_default_parameters
    @params['paired'] = dataset_has_column?('Read2')
  end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  end
  def next_dataset
    {'Name'=>@dataset['Name'], 
     'Report [Link]'=>File.join(@result_dir, "#{@dataset['Name']}/fastqc_report.html"), 
     'ReportDir [File]'=>File.join(@result_dir, "#{@dataset['Name']}")
    }
  end
  def commands
        command =<<-EOS
CORES=#{@params['cores']}
READS=#{@gstore_dir}/#{@dataset['Read1']}
NAME=#{@dataset['Name']}
FILENAME=${READS##*/}
FILEBASE=${FILENAME%%.*}
fastqc --extract -o . -t $CORES $READS
### fastqc generates the directory
FASTQC_DIR=${FILEBASE}_fastqc
mv $FASTQC_DIR $NAME
EOS
  end
end

