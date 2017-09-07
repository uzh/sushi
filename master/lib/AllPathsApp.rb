#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class AllPathsApp < SushiFabric::SushiApp
  def initialize
    super
@description =<<-EOS
Genome <i>deNovo</i> assembler based on  <a href="http://www.broadinstitute.org/software/allpaths-lg/blog/">AllPaths</a>.
The App requires a tab-separated input file with exactly 13 columns, i.e., the standard in_libs.csv columns (but still tab-sep, not comma-sep, as the original in_libs.csv would request) that
AllPaths requires with an additional column at the beginning called 'file'.
All the details about the in_libs.csv file and the various AllPaths options are described in the  <a href="ftp://ftp.broadinstitute.org/pub/crd/ALLPATHS/Release-LG/AllPaths-LG_Manual.pdf">AllPaths Manual</a>.
EOS
    @params['process_mode'] = 'DATASET'
    @name = 'AllPaths'
    @analysis_category = 'Assemble'
    @required_columns = ['file','library_name','project_name','organism_name','type','paired','frag_size','frag_stddev','insert_size','insert_stddev','read_orientation','genomic_start','genomic_end']
    @required_params = ['Estimated_Genome_Size','Estimated_Coverage_From_Fragment_Libraries','Estimated_Coverage_From_Jump_Libraries']
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
    @params['Ploidy','description'] = 'Is the genome to assemble haploid or diploid?' 
    @params['Remove_dodgy_reads'] = ['False','True']
    @params['Remove_dodgy_reads','description'] = 'Should AllPaths check for duplicates and low-quality reads? This is a rather quick-fix, a proper pre-processing step is recommended.'
    @params['Haplodify'] = ['False','True']
    @params['Haplodify', 'description'] = 'If the genome is diploid, activate this in case of an expected high level of heterozygosity.'
    @params['Additional_Options'] = ''
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def next_dataset
    {'Name'=>@dataset['Name'], 
     'Fasta [File]'=>File.join(@result_dir, "final.assembly.fasta"),
     'Fasta [File]'=>File.join(@result_dir, "final.contigs.fasta"),
     'Report [File]'=>File.join(@result_dir, "assembly.report"),
     'Stats [File]'=>File.join(@result_dir, "assembly_stats.report"),
    }.merge(extract_columns(@inherit_tags))
  end
#  def set_default_parameters
#    @params['refBuild'] = @dataset[0]['refBuild']
#  end

  def commands
    command =<<-EOS
ALL_PATHS_DIR=#{ALL_PATHS}
SIZE="#{@params['Estimated_Genome_Size']}"
F_COV="#{@params['Estimated_Coverage_From_Fragment_Libraries']}"
J_COV="#{@params['Estimated_Coverage_From_Jump_Libraries']}"
PLOIDY="#{@params['Ploidy']}"
CORES="#{@params['cores']}"
REM_DODGY="#{@params['Remove_dodgy_reads']}"
HAPLODIFY="#{@params['Haplodify']}"
TSV_FILE='#{@input_dataset_tsv_path}'
echo "group_name,library_name,file_name" > in_groups.csv
echo "library_name,project_name,organism_name,type,paired,frag_size,frag_stddev,insert_size,insert_stddev,read_orientation,genomic_start,genomic_end" > in_libs.csv
R --vanilla << EOT
x = read.table("$TSV_FILE", sep="\t", header=TRUE,blank.lines.skip = FALSE)
x1=subset(x,select=library_name:genomic_end)
write.table(x1, file="in_libs.csv", append=TRUE,sep=",", col.names = FALSE, row.names = FALSE, quote = FALSE) 

x2=cbind(subset(x,select=library_name),subset(x,select=library_name),subset(x,select=file))
write.table(x2, file="in_groups.csv",sep=",", append=TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
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
echo "$PLOIDY" > $(pwd)/AP_REF/AP_SUBDIR/ploidy

$ALL_PATHS_DIR/PrepareAllPathsInputs.pl DATA_DIR=$(pwd)/AP_REF/AP_SUBDIR  GENOME_SIZE=$SIZE FRAG_COVERAGE=$F_COV JUMP_COVERAGE=$J_COV  PLOIDY=$PLOIDY HOSTS=$CORES
$ALL_PATHS_DIR/RunAllPathsLG  PRE=$(pwd) REFERENCE_NAME=AP_REF RUN=ATT1 DATA_SUBDIR=AP_SUBDIR  HAPLOIDIFY=$HAPLODIFY TARGETS=standard  REMOVE_DODGY_READS_JUMP=$REM_DODGY
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

