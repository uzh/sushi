class JobMonitoringController < ApplicationController
  def index
    @job_list = `public/wfm_job_list`
    @job_list = @job_list.split(/\n/).map{|job| job.split(/,/)}
    @results  = {}
    @job_list.each do |job|
      if result = JobLog.find_by_job_id(job[0].to_i) and link=result.result_link.to_s and job[1] =~ /success/
        @results[job[0]]=link
      else
        @results[job[0]]=''
      end
    end
  end
  def print_log
    text = `public/wfm_getlog #{params[:job_id]} :with_err`
    render :text => text.gsub(/\n/,'<br />')
  end
end
