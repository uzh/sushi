#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class VariantCallerApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'VariantCaller'
    @analysis_category = 'Variant_Analysis'
@description =<<-EOS
Variant caller and variant annotator starting from a bam file. 
For calling variants, one can choose to using <a href="http://samtools.sourceforge.net/samtools.shtml">samtools+mpileup+bcftools</a> or <a href="http://www.broadinstitute.org/gatk/">GATK</a>.
To annotate variants, <a href="http://snpeff.sourceforge.net">snpEFF</a> is used. Please check <a href="http://snpeff.sourceforge.net/download.html#databases">here</a> whether the desired snpEFF database is avilable and needs to be downloaded.   
EOS
    @required_columns = ['Name','BAM','BAI', 'build']
    @required_params = ['min_depth_to_call_variants','build']
    # optional params
    @params['cores'] = '4'
    @params['ram'] = '10'
    @params['scratch'] = '100'
    @params['build'] = ref_selector
    @params['build','description'] = 'If human, then ensure that the build is hg_19_karyotypic'
    @params['snpEff_annotation'] = ['true','false']
    @params['snpEff_annotation','description'] = 'Annotate the variants? If yes, choose a snpEff database.'
#    @params['snpEff_database'] = {'select'=>''} 
#    @params['snpEff_database','description'] = 'If the database is not listed,  please check whether it is avilable and needs to be downloaded (link above).' 
#    Dir["/usr/local/ngseq/src/snpEff_v3.4/data/*"].sort.select{|build| File.directory?(build)}.each do |dir|
#      @params['snpEff_database'][File.basename(dir)] = File.basename(dir)
#    end
    @params['snpCaller'] = ['mpileup_bcftools','gatk']
    @params['snpCaller','description'] = 'Choose bewteen samtools+mpileup+bcftools and GATK. GATK is particularly recommended for human samples and cohort studies.'
    @params['paired'] = ['true','false']
    @params['paired','description'] = 'Are the reads paired?'
    @params['mpileupOptions'] = ''
    @params['bcftoolsOptions'] = ''
    @params['gatk_glm'] = ['SNP','INDEL','BOTH']
    @params['gatkOptions'] = '-baqGOP 30 -minIndelCnt 8 --min_base_quality_score 15 -stand_call_conf 15'
  end
  def next_dataset
    {'Name'=>@dataset['Name'], 
     'VCF [File]'=>File.join(@result_dir, "#{@dataset['Name']}.vcf"),
     'Gene_summary [File]'=>File.join(@result_dir, "#{@dataset['Name']}.genes.txt"),
     'Html [Link,File]'=>File.join(@result_dir, "#{@dataset['Name']}.html"),
     'build'=>@params['build']
    }.merge factor_dataset
  end
  def set_default_parameters
    @params['build'] = @dataset[0]['build']
    @params['min_depth_to_call_variants'] = '19'
  end

  def commands
    command =<<-EOS
set -e
set -o pipefail 
SAMTOOLS=#{GlobalVariables::SAMTOOLS}
BCFTOOLS=#{GlobalVariables::BCFTOOLS}
GATK_DIR=#{GlobalVariables::GATK_DIR}
PICARD_DIR=#{GlobalVariables::PICARD_DIR}
SNPEFF_DIR=#{GlobalVariables::SNPEFF_DIR}
SNP_CALLER="#{@params['snpCaller']}"
#SNPEFF_DATABASE="#{@params['snpEff_database']}"
MPILEUP_OPTIONS="#{@params['mpileupOptions']}"
BCF_OPTIONS="#{@params['bcftoolsOptions']}"
CORES="#{@params['cores']}"
GATK_GLM="#{@params['gatk_glm']}"
GATK_OPTIONS="#{@params['gatkOptions']}"
HSD=#{GlobalVariables::HUMAN_SNP_DATABASES}
MIN_DEPTH="#{@params['min_depth_to_call_variants']}"
PAIRED="#{@params['paired']}"
ANN="#{@params['snpEff_annotation']}"
BUILD="#{@params['build']}"

REF=/srv/GT/reference/#{@params['build']}/../../Sequence/WholeGenomeFasta/genome
MY_BAM=internal_grouped.lex.bam

if [ $PAIRED == "true" ]; then 
$SAMTOOLS view -F 4 -hb #{File.join(@gstore_dir, @dataset['BAM'])} | $SAMTOOLS rmdup - internal.nodup.bam
$SAMTOOLS index internal.nodup.bam
else
$SAMTOOLS view -F 4 -hb #{File.join(@gstore_dir, @dataset['BAM'])} | $SAMTOOLS rmdup -s - internal.nodup.bam
$SAMTOOLS index internal.nodup.bam 
fi

### SORT OUT GROUPS ISSUES ###
java -jar $PICARD_DIR/AddOrReplaceReadGroups.jar I=internal.nodup.bam \
   O=internal_grouped.lex.bam SORT_ORDER=coordinate RGID=ID_NAME TMP_DIR=/scratch \
   RGLB=Paired_end RGPL=illumina RGSM=project RGPU=BIOSEQUENCER
   MY_BAM=internal_grouped.lex.bam
$SAMTOOLS index $MY_BAM
   

### USE MPILEUP ####
 if [ $SNP_CALLER == 'mpileup_bcftools' ]; then  
  ### DETECTING VARIANTS BCFTOOLS ###
  $SAMTOOLS mpileup $MPILEUP_OPTIONS -uf $REF.fa $MY_BAM | $BCFTOOLS view -bvcg - > internal.bcf  
  $BCFTOOLS view  $BCF_OPTIONS  internal.bcf  > final.output.vcf 
 else

 ### USE GATK ###
 if [ ! -f $REF.fa.fai ]; then 
 $SAMTOOLS faidx $REF.fa
 fi  
 if [ ! -f $REF.dict ]; then
 java -jar $PICARD_DIR/CreateSequenceDictionary.jar R=$REF.fa O=$REF.dict
 fi 

#     if [ $REF == "/srv/GT/reference/#{@params['build']}/../../Sequence/WholeGenomeFasta/genome" ]; then
      human=$(grep "Homo_S" "/srv/GT/reference/#{@params['build']}/../../Sequence/WholeGenomeFasta/genome")
      if [ -n "$human"  ]; then
     ### HUMAN BEST PRACTICES ####

     ########FINDING POSSIBLE INDELS ####
     java -Xmx2g -jar $GATK_DIR/GenomeAnalysisTK.jar -T RealignerTargetCreator \
     -R $REF.fa -o paired_end.intervals --num_threads 4 -known $HSD/dbsnp_138.hg19.2.vcf \
     -I $MY_BAM

     ### REALINGING AROUND POSSIBLE INDELS ###
     java -Xmx4g -jar $GATK_DIR/GenomeAnalysisTK.jar   -I $MY_BAM  \
     -R $REF.fa -T IndelRealigner -known $HSD/dbsnp_138.hg19.2.vcf  \
     -targetIntervals paired_end.intervals -o $MY_BAM.real.trans.bam

     ### BASE RECALIBRATION ###
     java -Xmx4g -jar $GATK_DIR/GenomeAnalysisTK.jar \
     -T BaseRecalibrator \
     -I $MY_BAM.real.trans.bam \
     -R $REF.fa \
     -knownSites $HSD/dbsnp_138.hg19.2.vcf \
     -o recal_data.table

     ### APPLY BASE RECALIBRATION ###
     java -Xmx4g -jar $GATK_DIR/GenomeAnalysisTK.jar \
     -T PrintReads \
     -R $REF.fa \
     -I $MY_BAM.real.trans.bam \
     -BQSR recal_data.table \
     -o $MY_BAM.real.bam

     ### DETECTING VARIANTS GATK  ###
     java -jar $GATK_DIR/GenomeAnalysisTK.jar \
     -I $MY_BAM.real.bam  -log gatk_log.txt -nt $CORES \
     -o internal.vcf -R $REF.fa -T UnifiedGenotyper \
     -glm $GATK_GLM $GATK_OPTIONS --dbsnp  $HSD/dbsnp_138.hg19.2.vcf

     ### FILTER VARIANTS ###
     java -jar $GATK_DIR/GenomeAnalysisTK.jar \
     -T VariantFiltration \
     -R $REF.fa \
     -V internal.vcf \
     --filterExpression "DP < $MIN_DEPTH" \
     --filterName "my_snp_filter" \
     -o internal2.vcf

     grep "#" internal2.vcf  > internal3.vcf 
     grep "PASS" internal2.vcf   >> internal3.vcf


     ### RECALIBRATE VARIANTS ####
     java -Xmx16g -jar $GATK_DIR/GenomeAnalysisTK.jar \
     -T VariantRecalibrator \
     -R $REF.fa \
     -input  internal3.vcf \
     -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HSD/hapmap_3.3.hg19.reord.vcf \
     -resource:omni,known=false,training=true,truth=false,prior=12.0 $HSD/1000G_omni2.5.hg19.vcf \
     -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $HSD/dbsnp_138.hg19.2.vcf \
     -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
     -mode $GATK_GLM \
     -recalFile output.recal --num_threads 4 \
     -tranchesFile output.tranches \
     -rscriptFile output.plots.R -rf BadCigar 
 
     ### APPLY RECALIBRATION #####
     java -Xmx16g -jar $GATK_DIR/GenomeAnalysisTK.jar \
     -T ApplyRecalibration \
     -R $REF.fa \
     -input internal3.vcf \
     --ts_filter_level 99.0 \
     -tranchesFile output.tranches \
     -recalFile output.recal \
     -mode  $GATK_GLM --num_threads 4 \
     -o final.output.vcf -rf BadCigar 

     else 

     ### GENOME OTHER THAN HUMAN ###
     ########FINDING POSSIBLE INDELS ####
     java -Xmx2g -jar $GATK_DIR/GenomeAnalysisTK.jar -T RealignerTargetCreator \
     -R $REF.fa -o paired_end.intervals --num_threads 4 \
     -I $MY_BAM

     ### REALINGING AROUND POSSIBLE INDELS ###
     java -Xmx4g -jar $GATK_DIR/GenomeAnalysisTK.jar   -I $MY_BAM  \
     -R $REF.fa -T IndelRealigner \
     -targetIntervals paired_end.intervals -o $MY_BAM.real.bam

     ### DETECTING VARIANTS GATK  ###
     java -jar $GATK_DIR/GenomeAnalysisTK.jar \
     -I $MY_BAM.real.bam  -log gatk_log.txt -nt $CORES \
     -o internal.vcf -R $REF.fa -T UnifiedGenotyper \
     -glm $GATK_GLM $GATK_OPTIONS

     ### FILTER VARIANTS ###
     java -jar $GATK_DIR/GenomeAnalysisTK.jar \
     -T VariantFiltration \
     -R $REF.fa \
     -V internal.vcf \
     --filterExpression "DP < $MIN_DEPTH" \
     --filterName "my_snp_filter" \
     -o internal2.vcf

     grep "#" internal2.vcf  >  final.output.vcf
     grep "PASS" internal2.vcf   >> final.output.vcf
fi
fi

### ANNOTATION ####
if [ $ANN == "true" ]; then 
snpEffDir=/srv/GT/reference/"#{@params['build']}"/Genes/snpEff
mkdir -p $snpEffDir
awk -v str="DIRECTORY_FOR_DATA" -v str2="/srv/GT/reference/#{@params['build']}/Genes/snpEff" '{sub(str,str2,$0); print }' > $snpEffDir/snpEff.config

   ## CHECK IF DATABASE EXISTS AND CREATE IT IF NOT ###
    while [ ! -f "$snpEffDir/*/snpEffectPredictor.bin"  ]
     do
     if
     [ -f "$snpEffDir/temp.txt" ]; then
     sleep 1m
     continue
     else
     echo "database under contrsuction" > $snpEffDir/temp.txt
     base=$(echo /srv/GT/reference/"#{@params['build']}" | sed s/"\/"/" "/g | awk '{print $1}')
     provider=$(echo /srv/GT/reference/"#{@params['build']}" | sed s/"\/"/" "/g | awk '{print $2}')
     mkdir $snpEffDir
     mkdir $snpEffDir/$base.$provider
     echo "# $base"    >> $snpEffDir/snpEff.config
     echo "$provider.genome : $base"    >> $snpEffDir/snpEff.config
     echo "$provider.reference : $REF.fa"    >> $snpEffDir/snpEff.config
     ln $REF.fa $snpEffDir/$base.$provider/sequences.fa
     ln /srv/GT/reference/"#{@params['build']}"/Genes/genes.gtf $snpEffDir/sequences.fa/$base.$provider
     java -Xmx2g -jar $SNPEFF_DIR/snpEff.jar build -c $snpEffDir/snpEff.config -gtf22 -v "$base.$provider"
    fi
   done

    java -Xmx2g -jar $SNPEFF_DIR/snpEff.jar -s #{@dataset['Name']}.html -c $snpEffDir/snpEff.config $base.$provider -v final.output.vcf  > #{@dataset['Name']}.vcf 
fi
EOS
    command
  end
end

if __FILE__ == $0

end
