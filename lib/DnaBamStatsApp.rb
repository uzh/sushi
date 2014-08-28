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
    @required_columns = ['Name','BAM']
    @required_params = ['cores','ram','build','sortedBam']
    @params['cores'] = '8'
    @params['ram'] = '60'
    @params['scratch'] = '100'
    @params['build'] = ref_selector
    @params['sortedBam'] = ['false','true']
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
    @params['build'] = @dataset[0]['build']
  end
  def commands
    command =<<-EOS
REF=/srv/GT/reference/#{@params['build']}/../../Sequence/WholeGenomeFasta/genome.fa
SAMTOOLS=#{GlobalVariables::SAMTOOLS}
SAMSTAT=#{GlobalVariables::SAMSTAT}
QUALIMAP=#{GlobalVariables::QUALIMAP}
PICARD_DIR=#{GlobalVariables::PICARD_DIR}
###samstat
$SAMTOOLS view -ub #{File.join(@gstore_dir, @dataset['BAM'])} | $SAMSTAT -f bam -n #{@dataset['Name']}.samstat
###qualimap
unset DISPLAY ; $QUALIMAP bamqc -bam  #{File.join(@gstore_dir, @dataset['BAM'])} -c -nt #{@params['cores']} --java-mem-size=10G -outdir #{@dataset['Name']}
rm #{@dataset['Name']}/coverage.txt
###picard 
/usr/bin/java -Xmx10g -jar $PICARD_DIR/CollectGcBiasMetrics.jar OUTPUT=#{@dataset['Name']}.gc.dat SUMMARY_OUTPUT=#{@dataset['Name']}.gc.sum.dat INPUT=#{File.join(@gstore_dir, @dataset['BAM'])} CHART_OUTPUT=#{@dataset['Name']}.picard.pdf ASSUME_SORTED=#{@params['sortedBam']} REFERENCE_SEQUENCE=$REF VALIDATION_STRINGENCY=LENIENT
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

