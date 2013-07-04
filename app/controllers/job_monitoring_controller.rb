class JobMonitoringController < ApplicationController
  def index
    @job_list = `public/wfm_job_list -d #{WORKFLOW_MANAGER} -p #{session[:project]}`
    @job_list = @job_list.split(/\n/).map{|job| job.split(/,/)}
    @total = @job_list.length

    # pager
    @page_unit = if page = params[:page] and unit = page[:unit]
                   session[:job_page_unit] = unit.to_i
                 elsif unit = session[:job_page_unit]
                   unit.to_i
                 else
                   session[:job_page_unit] = 10
                 end
    current_page = params[:format]
    @current_page = (current_page||1).to_i
    @page_list = (1..(@job_list.length.to_f/@page_unit).ceil).to_a
    start = (@current_page - 1) * @page_unit
    last  = @current_page * @page_unit - 1
    @job_list = @job_list[start..last]
    @submit_jobs = []
    @job_list.each_with_index do |job, i|
      if submit_job = Job.find_by_submit_job_id(job[0].to_i)
        @submit_jobs[i] = submit_job
      end
    end
  end
  def print_log
    text = `public/wfm_get_log #{params[:job_id]} :with_err #{WORKFLOW_MANAGER}`
    render :text => text.gsub(/\n/,'<br />')
  end
  def print_script
    text = `public/wfm_get_script #{params[:job_id]} #{WORKFLOW_MANAGER}`
    render :text => text.gsub(/\n/,'<br />')
  end
end
