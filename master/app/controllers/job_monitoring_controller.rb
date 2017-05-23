class JobMonitoringController < ApplicationController
  def index
    public_dir = File.expand_path('../../../public', __FILE__)
    @page_unit = 100
    @job_list = if option=params[:option] and option[:all_job_list] 
                  @page_unit = 1000
                  @all_job_list=true
                  session[:all_job_list] = true
                  #@@workflow_manager.job_list(false, nil)
                  command = "wfm_job_list -d #{SushiFabric::WORKFLOW_MANAGER}"
                  `#{command}`
                elsif option=params[:option] and option[:project_job_list] 
                  @all_job_list=false
                  session[:all_job_list] = false
                  #@@workflow_manager.job_list(false, session[:project])
                  command = "wfm_job_list -p #{session[:project]} -d #{SushiFabric::WORKFLOW_MANAGER}"
                  `#{command}`
                elsif session[:all_job_list]
                  @page_unit = 1000
                  @all_job_list=true
                  session[:all_job_list] = true
                  command = "wfm_job_list -d #{SushiFabric::WORKFLOW_MANAGER}"
                  `#{command}`
                else
                  @all_job_list=false
                  session[:all_job_list] = false
                  command = "wfm_job_list -p #{session[:project]} -d #{SushiFabric::WORKFLOW_MANAGER}"
                  `#{command}`
                end
    @job_list = @job_list.split(/\n/).map{|job| job.split(/,/)}
    @total = @job_list.length

    # pager
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
    text = @@workflow_manager.get_log(params[:job_id], :with_err)
    render :text => text.gsub(/\n/,'<br />')
  end
  def print_script
    text = 'no script found'
    if sushi_job_id = params[:sushi_job_id] and
      job = Job.find_by_id(sushi_job_id.to_i) and
      script_path = job.script_path and File.exist?(script_path)
      text = File.read(script_path)
    else
      text = @@workflow_manager.get_script(params[:job_id])
    end
    render :text => text.gsub(/\n/,'<br />')
  end
  def kill_job
    @status = 'kill job failed'
    if @job_id = params[:id]
      public_dir = File.expand_path('../../../public', __FILE__)
      @status = @@workflow_manager.kill_job(@job_id)
      @command = "wfm_kill_job -i #{@job_id} -d #{SushiFabric::WORKFLOW_MANAGER}"
    end
  end
  def change_status
    if @job_id = params[:id]
      public_dir = File.expand_path('../../../public', __FILE__)
      status = @@workflow_manager.status(@job_id)
      if status and @status = status.split(',').first
        if @status == 'success'
          @@workflow_manager.status(@job_id, "fail")
        elsif @status == 'fail'
          @@workflow_manager.status(@job_id, "success")
        end
      end
    end
    redirect_to :controller => "job_monitoring"
  end
end
