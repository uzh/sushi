#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class GatkRnaApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'GatkRna'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Variants'
    @description =<<-EOS
Variant calling for RNA-seq using <a href='https://software.broadinstitute.org/gatk/'>GATK Queue</a> and recommended filtering.
    EOS
    @required_columns = ['Name', 'BAM', 'BAI', 'refBuild', 'Species']
    @required_params = ['name', 'paired']
    @params['threads'] = '16'
    @params['threads', 'description'] = 'number of data threads for jobs that support -nt'
    @params['cthreads'] = '16'
    @params['cthreads', 'description'] = 'number of CPU threads for jobs that support -nct; scatter jobs will get cthreads/sj threads each'
    @params['sj'] = '16'
    @params['sj', 'description'] = 'number of scatter jobs to generate'
    @params['maxConcurrentRun'] = '20'
    @params['maxConcurrentRun', 'description'] = 'maximum number of jobs to schedule concurrently'
    @params['ram'] = '100'
    @params['scratch'] = '200'
    @params['paired'] = false
    @params['name'] = 'GATK_RnaVariants'
    @params['refBuild'] = ref_selector
    @params['rgpl'] = 'illumina'
    @params['rgpl', 'description'] = 'platform/technology used to produce reads (used to define read groups)'
	 @params['dontUseSoftClippedBases'] = ['true', 'false']
    @params['mail'] = ""
    @modules = ["Tools/samtools", "Dev/jdk", "Variants/Queue", "Tools/Picard", "Dev/R"]
  end
  def next_dataset
    {
     'Name'=>@params['name'],
     'VCF [File]'=>File.join(@result_dir, "#{@params['name']}.vcf.gz"),
     'VCF_index [File]'=>File.join(@result_dir, "#{@params['name']}.vcf.gz.tbi"),
     'Species'=>(dataset = @dataset.first and dataset['Species']),
     'refBuild'=>@params['refBuild']
    }
  end
  def set_default_parameters
    @params['cores'] = '1'
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('paired')
      @params['paired'] = @dataset[0]['paired']
    end
  end

  def commands
    run_RApp('EzAppGatkRna')
  end
end

if __FILE__ == $0

end
