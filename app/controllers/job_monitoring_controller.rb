class JobMonitoringController < ApplicationController
  def index
    public_dir = File.expand_path('../../../public', __FILE__)
    @job_list = if option=params[:option] and option[:all_job_list]
                  @all_job_list=true
                  `#{public_dir}/wfm_job_list -d #{SushiFabric::WORKFLOW_MANAGER}`
                else
                  `#{public_dir}/wfm_job_list -d #{SushiFabric::WORKFLOW_MANAGER} -p #{session[:project]}`
                end
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
    public_dir = File.expand_path('../../../public', __FILE__)
    text = `#{public_dir}/wfm_get_log #{params[:job_id]} :with_err #{SushiFabric::WORKFLOW_MANAGER}`
    render :text => text.gsub(/\n/,'<br />')
  end
  def print_script
    public_dir = File.expand_path('../../../public', __FILE__)
    text = `#{public_dir}/wfm_get_script #{params[:job_id]} #{SushiFabric::WORKFLOW_MANAGER}`
    render :text => text.gsub(/\n/,'<br />')
  end
  def kill_job
    @status = 'kill job failed'
    if @job_id = params[:id]
      public_dir = File.expand_path('../../../public', __FILE__)
      @command = "#{public_dir}/wfm_kill_job -i #{@job_id} -d #{SushiFabric::WORKFLOW_MANAGER}"
      @status = `#{@command}`
      @command = "wfm_kill_job -i #{@job_id} -d #{SushiFabric::WORKFLOW_MANAGER}"
    end
  end
  def change_status
    if @job_id = params[:id]
      public_dir = File.expand_path('../../../public', __FILE__)
      @command = "#{public_dir}/wfm_status #{@job_id} #{SushiFabric::WORKFLOW_MANAGER}"
      @status = `#{@command}`.split(',').first
      if @status == 'success'
        @command = "#{public_dir}/wfm_status #{@job_id} #{SushiFabric::WORKFLOW_MANAGER} fail"
        `#{@command}`
      elsif @status == 'fail'
        @command = "#{public_dir}/wfm_status #{@job_id} #{SushiFabric::WORKFLOW_MANAGER} success"
        `#{@command}`
      end
    end
    redirect_to :controller => "job_monitoring"
  end
end
