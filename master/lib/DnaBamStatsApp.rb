#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class DnaBamStatsApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'DNA BamStats'
    @analysis_category = 'QC'
    @description =<<-EOS
Runs the following tools that check alignment statistics for the DNA applications<br/>
<a href='http://samstat.sourceforge.net/'>SAMStat</a><br/>
<a href='http://qualimap.bioinfo.cipf.es/'>QualiMap</a><br/>
<a href='http://broadinstitute.github.io/picard/command-line-overview.html#CollectGcBiasMetrics'>picard/CollectGcBiasMetrics</a><br/>
EOS
    @required_columns = ['Name','BAM']
    @required_params = ['cores','ram','refBuild','sortedBam']
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '30'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '100'
    @params['scratch', "context"] = "slurm"
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "referfence genome assembly"
    @params['sortedBam'] = true
    @params['sortedBam', "context"] = "DnaBamStats"
    @modules = ["Dev/R", "Tools/Picard", "QC/Qualimap", "QC/SAMStat", "Tools/samtools", "Dev/jdk/8"]
  end
  def next_dataset
    samstat_link = File.join(@result_dir, "#{@dataset['Name']}.samstat.html")
    qualimap_link = File.join(@result_dir, "#{@dataset['Name']}", 'qualimapReport.html')
    picard_link = File.join(@result_dir, "#{@dataset['Name']}.picard.pdf")
    {'Name'=>@dataset['Name'],
     'Samstat Result [File]'=>samstat_link,
     'Qualimap Result [File]'=>File.join(@result_dir, "#{@dataset['Name']}"),
     'Picard Result [File]'=>picard_link,
     'Samstat Report [Link]'=>samstat_link,
     'Qualimap Report [Link]'=>qualimap_link,
     'Picard Report  [Link]'=>picard_link
    }
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
  end
  def commands
    command =<<-EOS
#echo Display=$DISPLAY
set -x
REF=/srv/GT/reference/#{@params['refBuild']}/../../Sequence/WholeGenomeFasta/genome.fa
BAM_FILE=#{File.join(@gstore_dir, @dataset['BAM'])}
FN=$(basename $BAM_FILE)
###samstat
cp -a $BAM_FILE .
samstat $FN
mv {$FN,#{@dataset['Name']}}.samstat.html
###qualimap
unset DISPLAY ; qualimap bamqc -bam $FN -c -nt #{@params['cores']} --java-mem-size=#{@params['ram']}G -outdir #{@dataset['Name']} > qualimap.out
###picard 
/usr/local/ngseq/packages/Dev/jdk/21/bin/java -Xmx#{@params['ram']}g -jar $Picard_jar CollectGcBiasMetrics OUTPUT=#{@dataset['Name']}.gc.dat SUMMARY_OUTPUT=#{@dataset['Name']}.gc.sum.dat INPUT=$FN CHART_OUTPUT=#{@dataset['Name']}.picard.pdf ASSUME_SORTED=#{@params['sortedBam']} REFERENCE_SEQUENCE=$REF VALIDATION_STRINGENCY=LENIENT > picard.out 2> picard.err
EOS
    command
  end
end

if __FILE__ == $0
  usecase = DnaBamStats.new

  usecase.project = "p1001"
  usecase.user = 'qiwei'

  # set user parameter
  # for GUI sushi
  #usecase.params['process_mode'].value = 'SAMPLE'
  usecase.params['cores'] = 1
  usecase.params['node'] = 'fgcz-c-048'

  # also possible to load a input dataset from Sushi DB
  usecase.dataset_sushi_id = 1

  # run (submit to workflow_manager)
  #usecase.run
  usecase.test_run

end

