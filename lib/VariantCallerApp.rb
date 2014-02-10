#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'

class VariantCallerApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'VariantCaller'
    @analysis_category = 'Variant_Analysis'
    @required_columns = ['Name','BAM','BAI', 'Build']
		@required_params = []
    # optional params
    @params['cores'] = '4'
    @params['ram'] = '10'
    @params['scratch'] = '100'
    @params['glm'] = ['BOTH', 'INDEL', 'SNP'] 
    @params['build'] = {'select'=>''}
    Dir["/srv/GT/reference/*/*/*"].sort.select{|build| File.directory?(build)}.each do |dir|
      @params['build'][dir.gsub(/\/srv\/GT\/reference\//,'')] = File.basename(dir)
    end
    @params['gatkOptions'] = '-minIndelCnt 8 --min_base_quality_score 15 -stand_call_conf 15 -baqGOP 30'
  end
  def next_dataset
    {'Name'=>@dataset['Name'], 
     'VCF [File]'=>File.join(@result_dir, "#{@dataset['Name']}.vcf"), 
     'Build'=>@dataset['Build']
    }
  end
  def commands
    command =<<-EOS
samtools rmdup #{File.join(@gstore_dir, @dataset['BAM'])} internal.nodup.bam
samtools index internal.nodup.bam

### SORT OUT GROUPS ISSUES ###
java -jar /usr/local/ngseq/stow/picard-tools-1.96/bin/AddOrReplaceReadGroups.jar I=internal.nodup.bam \
O=internal_grouped.bam SORT_ORDER=coordinate RGID=ID_NAME TMP_DIR=/scratch \
RGLB=Paired_end RGPL=illumina RGSM=project RGPU=BIOSEQUENCER

### REINDEXING CORRECTLY GROUPED FILES ###
samtools index internal_grouped.bam

### DETECTING VARIANTS ###
GATK_DIR=/usr/local/ngseq/src/GenomeAnalysisTK-2.8-1-g932cd3a
REF=$(ls -l /srv/GT/reference/*/*/*/Sequence/WholeGenomeFasta/genome.fa | #{@params['build']} | awk '{print $9}')

java -Xmx4g -jar $GATK_DIR/GenomeAnalysisTK.jar \
  -I internal_grouped.bam  -log gatk_log.txt -nt #{@params['cores']} \
  -o internal.vcf -R $REF -T UnifiedGenotyper \
  -glm #{@params['glm']} \
  #{@params['gatkOptions']}
### ANNOTATION ####
SNPEFF_DIR=/usr/local/ngseq/src/snpEff_v3.4/
#cd $SNPEFF_DIR
java -Xmx4g -jar $SNPEFF_DIR/snpEff.jar -c $SNPEFF_DIR/snpEff.config  athalianaTair10  -v internal.vcf > #{@dataset['Name']}.vcf
EOS
    command
  end
end

if __FILE__ == $0
  usecase = VariantCallerApp.new

  usecase.project = "p1001"
  usecase.user = 'masamasa'

  # set user parameter
  # for GUI sushi
  #usecase.params['process_mode'].value = 'SAMPLE'
  usecase.params['cores'] = 8
  usecase.params['node'] = 'fgcz-c-048'

  # also possible to load a parameterset csv file
  # mainly for CUI sushi
  #usecase.parameterset_tsv_file = 'tophat_parameterset.tsv'
  #usecase.parameterset_tsv_file = 'test.tsv'

  # set input dataset
  # mainly for CUI sushi
  #usecase.dataset_tsv_file = 'tophat_dataset.tsv'

  # also possible to load a input dataset from Sushi DB
  usecase.dataset_sushi_id = 1

  # run (submit to workflow_manager)
  #usecase.run
  usecase.test_run

end

