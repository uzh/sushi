class RunScriptController < ApplicationController
  def index
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
