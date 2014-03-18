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
    @required_columns = ['Name','BAM','BAI', 'build']
    @required_params = []
    # optional params
    @params['cores'] = '4'
    @params['ram'] = '10'
    @params['scratch'] = '100'
    @params['build'] = ref_selector
    @params['snpEff_database'] = {'select'=>''}
    Dir["/usr/local/ngseq/src/snpEff_v3.4/data/*"].sort.select{|build| File.directory?(build)}.each do |dir|
      @params['snpEff_database'][File.basename(dir)] = File.basename(dir)
    end
    @params['snpCaller'] = ['mpileup_bcftools','gatk']
    @params['mpileupOtions'] = ''
    @params['bcftoolsOtions'] = ''
    @params['gatk_glm'] = ['SNP','INDEL','BOTH']
    @params['gatkOptions'] = '-baqGOP 30 -minIndelCnt 8 --min_base_quality_score 15 -stand_call_conf 15'
  end
  def next_dataset
    {'Name'=>@dataset['Name'], 
     'VCF [File]'=>File.join(@result_dir, "#{@dataset['Name']}.vcf"),
     'Html [Link,File]'=>File.join(@result_dir, "#{@dataset['Name']}.html"),
     'build'=>@params['build']
    }
  end
  def set_default_parameters
    @params['build'] = @dataset[0]['build']
  end

  def commands
    command =<<-EOS
SAMTOOLS=#{GlobalVariables::SAMTOOLS}
BCFTOOLS=#{GlobalVariables::BCFTOOLS}
GATK_DIR=#{GlobalVariables::GATK_DIR}
PICARD_DIR=#{GlobalVariables::PICARD_DIR}
SNPEFF_DIR=#{GlobalVariables::SNPEFF_DIR}
SNP_CALLER=#{@params['snpCaller']}
SNPEFF_DATABASE=#{@params['snpEff_database']}
MPILEUP_OPTIONS=#{@params['mpileupOtions']}
BCF_OPTIONS=#{@params['bcftoolsOtions']}
CORES=#{@params['cores']}
GATK_GLM=#{@params['gatk_glm']}
GATK_OPTIONS=#{@params['gatkOptions']}

REF=/srv/GT/reference/#{@params['build']}/../../Sequence/WholeGenomeFasta/genome.fa
MY_BAM=internal_grouped.lex.bam

$SAMTOOLS view -F 4 -hb #{File.join(@gstore_dir, @dataset['BAM'])} | $SAMTOOLS rmdup - internal.nodup.bam
$SAMTOOLS index internal.nodup.bam

### SORT OUT GROUPS ISSUES ###
java -jar $PICARD_DIR/AddOrReplaceReadGroups.jar I=internal.nodup.bam \
   O=internal_grouped.lex.bam SORT_ORDER=coordinate RGID=ID_NAME TMP_DIR=/scratch \
   RGLB=Paired_end RGPL=illumina RGSM=project RGPU=BIOSEQUENCER
   BAMFILE=$MY_BAM 
$SAMTOOLS index $MY_BAM
   
if [ $SNP_CALLER == 'mpileup_bcftools' ]; then  
 ### DETECTING VARIANTS BCFTOOLS ###
 $SAMTOOLS mpileup $MPILEUP_OPTIONS -uf $REF $MY_BAM | $BCFTOOLS view -bvcg - > internal.bcf  
 $BCFTOOLS view  $BCF_OPTIONS  internal.bcf  > internal.vcf 
else
  ### DETECTING VARIANTS GATK  ###
  java -jar $GATK_DIR/GenomeAnalysisTK.jar \
   -I $MY_BAM  -log gatk_log.txt -nt $CORES \
  -o internal.vcf -R $REF -T UnifiedGenotyper \
  -glm $GATK_GLM $GATK_OPTIONS
fi
### ANNOTATION ####
java -Xmx2g -jar $SNPEFF_DIR/snpEff.jar -s #{@dataset['Name']}.html -c $SNPEFF_DIR/snpEff.config $SNPEFF_DATABASE -v internal.vcf > #{@dataset['Name']}.vcf
EOS
    command
  end
end

if __FILE__ == $0

end

