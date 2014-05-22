#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class AllPathsDeNovoApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'AllPaths'
    @analysis_category = 'DeNovoAssembler'
    @required_columns = ['library','file','project','organism','type','frag_size','frag_stddev','insert_size','insert_stddev','read_orientation']
    @required_params = []
    # optional params
    @params['cores'] = '48'
    @params['ram'] = '300'
    @params['scratch'] = '500'
    @params['Path to in_libs.csv file'] = ''
    @params['Path to in_groups.csv file'] = ''
    @params['Estimated_Genome_Size'] = ''
    @params['Estimated_Coverage_From_Fragment_Libraries'] = ''
    @params['Estimated_Coverage_From_Jump_Libraries'] = ''
    @params['Ploidy']  = ['1','2']
    @params['Remove_dodgy_reads?'] = ['True','False']
    @params['Haplodify?'] = ['True','False']
    @params['Additional_Options'] = ''
  end
  def next_dataset
    {'Name'=>@dataset['Name'], 
     'Fasta [File]'=>File.join(@result_dir, "#{@dataset['Name']}.fasta"),
    }.merge factor_dataset
  end
  def set_default_parameters
    @params['build'] = @dataset[0]['build']
  end

  def commands
    command =<<-EOS
ALL_PATHS_DIR=#{GlobalVariables::ALL_PATHS}
SIZE="#{@params['Estimated_Genome_Size']}"
F_COV="#{@params['Estimated_Coverage_From_Fragment_Libraries']}"
J_COV="#{@params['Estimated_Coverage_From_Jump_Libraries']}"
PLOIDY="#{@params['Ploidy']}"
CORES="#{@params['cores']}"
REM_DODGY="#{@params['Remove_dodgy_reads?']}"
echo "group_name,library_name,file_name" > in_groups.csv
echo  #{@dataset['Name']} 
echo "library_name,project_name,organism_name,type,paired,frag_size,frag_stddev,insert_size,insert_stddev,read_orientation,genomic_start,genomic_end" > in_libs.csv

$ALL_PATHS_DIR/PrepareAllPathsInputs.pl DATA_DIR=$  GENOME_SIZE=$SIZE FRAG_COVERAGE=$F_COV JUMP_COVERAGE=$J_COV  PLOIDY=$PLOIDY HOSTS=$CORES
$ALL_PATHS_DIR/RunAllPathsLG  PRE=$WD REFERENCE_NAME=ALL_PATHS RUN=ATT1 DATA_SUBDIR=3_PE_AND_ALL_PROTON  HAPLOIDIFY=True TARGETS=standard  REMOVE_DODGY_READS_JUMP=$REM_DODGY
EOS
    command
  end
end

if __FILE__ == $0

end

