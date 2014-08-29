#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class AllPathsDeNovoApp < SushiFabric::SushiApp
  def initialize
    super
    @params['process_mode'] = 'DATASET'
    @name = 'AllPaths'
    @analysis_category = 'DeNovoAssembler'
   @required_columns = ['library','file','project','organism','type','paired','frag_size','frag_stddev','insert_size','insert_stddev','read_orientation','genomic_start','genomic_end']
    @required_params = []
    # optional params
    @params['cores'] = '48'
    @params['ram'] = '300'
    @params['scratch'] = '500'
    @params['Estimated_Genome_Size'] = ''
    @params['Estimated_Coverage_From_Fragment_Libraries'] = ''
    @params['Estimated_Coverage_From_Fragment_Libraries','description'] = 'A coverage between 50x and 70x for at least one of the libraries is recommended.'
    @params['Estimated_Coverage_From_Jump_Libraries'] = ''
    @params['Estimated_Coverage_From_Jump_Libraries','description'] = 'A coverage between 40x and 50x for at least one of the libraries is recommended.'
    @params['Ploidy']  = ['1','2']
    @params['Remove_dodgy_reads'] = ['True','False']
    @params['Haplodify'] = ['True','False']
    @params['Additional_Options'] = ''
  end
  def next_dataset
    {'Name'=>@dataset['Name'], 
     'Fasta [File]'=>File.join(@result_dir, "final.assembly.fasta"),
     'Fasta [File]'=>File.join(@result_dir, "final.contigs.fasta"),
     'Report [File]'=>File.join(@result_dir, "assembly.report"),
     'Stats [File]'=>File.join(@result_dir, "assembly_stats.report"),
    }.merge factor_dataset
  end
#  def set_default_parameters
#    @params['build'] = @dataset[0]['build']
#  end

  def commands
    command =<<-EOS
ALL_PATHS_DIR=#{GlobalVariables::ALL_PATHS}
SIZE="#{@params['Estimated_Genome_Size']}"
F_COV="#{@params['Estimated_Coverage_From_Fragment_Libraries']}"
J_COV="#{@params['Estimated_Coverage_From_Jump_Libraries']}"
PLOIDY="#{@params['Ploidy']}"
CORES="#{@params['cores']}"
REM_DODGY="#{@params['Remove_dodgy_reads?']}"
HAPLODIFY="#{@params['Haplodify?']}"
TSV_FILE='#{@input_dataset_tsv_path}'
echo "group_name,library_name,file_name" > in_groups.csv
echo "library_name,project_name,organism_name,type,paired,frag_size,frag_stddev,insert_size,insert_stddev,read_orientation,genomic_start,genomic_end" > in_libs.csv
R --vanilla << EOT
x = read.table("$TSV_FILE", sep="\t", header=TRUE,blank.lines.skip = FALSE)
x1=cbind(subset(x,select=library),subset(x,select=project:genomic_end))
write.table(x1, file="in_libs.csv", append=TRUE,sep=",", col.names = FALSE, row.names = FALSE, quote = FALSE) 

x2=cbind(subset(x,select=library),subset(x,select=library:file),sep=",")
write.table(x2, file="in_groups.csv", append=TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
EOT
sed -i s/"NA"/""/g in_groups.csv 
sed -i s/"NA"/""/g in_libs.csv
#echo "group_name,library_name,file_name" > in_groups.csv
#echo "#{@dataset['library']},#{@dataset['library']},#{@dataset['file']}" >> in_groups.csv 

#echo "library_name,project_name,organism_name,type,paired,frag_size,frag_stddev,insert_size,insert_stddev,read_orientation,genomic_start,genomic_end" > in_libs.csv
#echo "#{@dataset['library']},#{@dataset['project']},#{@dataset['organism']},#{@dataset['type']},1,#{@dataset['frag_size']},#{@dataset['frag_stddev']},#{@dataset['insert_size']}",#{@dataset['insert_stddev']},\
##{@dataset['read_orientation']},," >> in_libs.csv
cp in_* /srv/GT/analysis/p1438/AppTest
mkdir $(pwd)/AP_REF
mkdir $(pwd)/AP_REF/AP_SUBDIR 
echo "$PLOIDY" > AP_REF/AP_SUBDIR/ploidy

$ALL_PATHS_DIR/PrepareAllPathsInputs.pl DATA_DIR=$(pwd)/AP_REF/AP_SUBDIR  GENOME_SIZE=$SIZE FRAG_COVERAGE=$F_COV JUMP_COVERAGE=$J_COV  PLOIDY=$PLOIDY HOSTS=$CORES
$ALL_PATHS_DIR/RunAllPathsLG  PRE=$(pwd) REFERENCE_NAME=$(pwd)/AP_REF RUN=ATT1 DATA_SUBDIR=AP_SUBDIR  HAPLOIDIFY=$HAPLODIFY TARGETS=standard  REMOVE_DODGY_READS_JUMP=$REM_DODGY
mv AP_REF/AP_SUBDIR/ATT1/ASSEMBLIES/test/final.assembly.fasta ./
mv AP_REF/AP_SUBDIR/ATT1/ASSEMBLIES/test/final.contigs.fasta  ./
mv AP_REF/AP_SUBDIR/ATT1/ASSEMBLIES/test/assembly.report ./
mv AP_REF/AP_SUBDIR/ATT1/ASSEMBLIES/test/assembly_stats.report ./
EOS
    command
  end
end

if __FILE__ == $0

end

