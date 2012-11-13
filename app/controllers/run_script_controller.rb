class RunScriptController < ApplicationController
  def index
  end
  def run_sample
    @result = `sh public/sample.sh`
    render "run_script/run_sample"
  end
  def run_sample2
    render :text => 'start run_sample2.sh'
  end
end
