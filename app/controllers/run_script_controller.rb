$:.push '/srv/SushiFabric/plugins/bfabric/lib/'

require 'bfabric'
require 'tempfile'

class RunScriptController < ApplicationController
  def index
    @ext_job_id = 3506
#    @status = Bfabric.get_external_job_status @ext_job_id

    # search job_scripts
    @job_scripts = Dir['public/*.sh'].sort.to_a
    @data_sets = DataSet.all 
  end
  def set_parameters
    @data_sets = DataSet.all 
    @params = params
    @job_script = params[:script][:path]
    data_set_id = params[:dataset][:id]
    @data_set = @data_sets[data_set_id.to_i-1]
    @parameters = []
    File.readlines(@job_script).each do |line|
      if line =~ /#PARAMETER/
        x = line.split
        parameter = x[1]
        @parameters << parameter
      end
    end
    @default_parameters = {}
    File.readlines(@job_script).each do |line|
      @parameters.each do |parameter|
        if line =~ /#{parameter}=(.+)/
          @default_parameters[parameter]=$1
        end
      end
    end
    render "run_script/set_parameters"
  end
  def run_application
    @data_sets = DataSet.all 
    @params = params
    data_set_id = params[:dataset][:id]
    @data_set = @data_sets[data_set_id.to_i-1]
    @sample_ids = params[:sample_id].map{|i| i.to_i}.sort
    inputs = @data_set.data_lists.select{|data_list| @sample_ids.include?(data_list.sample.id)}
    @inputs_string = inputs.map{|data_list| data_list.sample.path}
    @job_script = params[:job_script][:path]
    parameters = params[:parameter]
    parameters['INPUT']="'"+@inputs_string.join(' ')+"'" unless @inputs_string.empty?
    script = Tempfile.open([File.basename(@job_script).split(/\./).first+'-','.sh'], 'public')
    File.readlines(@job_script).each do |line|
      flag = true
      parameters.keys.each do |parameter|
        if line =~ /#{parameter}=/
          script.print parameter, '=', parameters[parameter], "\n"
          flag = false
          break
        end
      end
      script.print line if flag
    end
    script.close
    sleep 1
    @job_id = `public/wfm_monitoring public/#{File.basename(script.path)}`
    script.delete
    render "run_script/submit_job"
  end

	def confirm
    render "run_script/confirm"
	end
	def run_fastqc
    @result = `ruby public/wfm_monitoring public/fastqc.sh druby://fgcz-s-024.fgcz-net.unizh.ch:12345`
		sleep 1
    render "run_script/run_sample"
	end
  def run_sample
    #@result = `sh public/sample.sh`
    #render "run_script/run_sample"
    #render "run_script/fastqc_report"
    render "public/fastqc_report"
  end
  def run_sample2
    render :text => 'start run_sample2.sh'
  end
end
