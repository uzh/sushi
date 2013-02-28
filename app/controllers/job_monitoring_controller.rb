class JobMonitoringController < ApplicationController
  def index
    @job_list = `public/wfm_job_list`
    @job_list = @job_list.split(/\n/).map{|job| job.split(/,/)}
    @results  = {}
    @job_list.each do |job|
      if job[1] =~ /success/ and result_link = `public/wfm_get_result_link #{job[0]}` 
        @results[job[0]]=result_link 
      else
        @results[job[0]]=''
      end
    end
  end
  def print_log
    text = `public/wfm_get_log #{params[:job_id]} :with_err`
    render :text => text.gsub(/\n/,'<br />')
  end
  def print_script
    text = `public/wfm_get_script #{params[:job_id]}`
    render :text => text.gsub(/\n/,'<br />')
  end
end
