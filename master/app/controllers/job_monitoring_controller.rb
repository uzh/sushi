class JobMonitoringController < ApplicationController
  def fetch_jobs(params, project_number)
    option = params[:option]

    if option && option[:all_job_list] # all projects
      configure_job_list(all_job_list: true, page_unit: 1000)
      Job.joins(data_set: :project)
         .where(projects: { number: project_number })
         .order(id: :desc)
         .limit(100)
    elsif option && option[:project_job_list] # project specific
      configure_job_list(all_job_list: false)
      Job.joins(data_set: :project)
         .order(id: :desc)
         .limit(100)
    elsif session[:all_job_list] # all projects from session
      configure_job_list(all_job_list: true, page_unit: 1000)
      Job.joins(data_set: :project)
         .where(projects: { number: project_number })
         .order(id: :desc)
         .limit(100)
    else # project specific
      configure_job_list(all_job_list: false)
      Job.joins(data_set: :project)
         .where(projects: { number: project_number })
         .order(id: :desc)
    end
  end
  def configure_job_list(all_job_list:, page_unit: nil)
    @all_job_list = all_job_list
    session[:all_job_list] = all_job_list
    @page_unit = page_unit if page_unit
  end

  def index
    public_dir = File.expand_path('../../../public', __FILE__)
    @page_unit = 100
    project_number = session[:project]

    jobs = fetch_jobs(params, project_number)

    @job_list = jobs.map{|job| [job.id, job.status.to_s, "job_script", job.start_time ? "#{job.start_time.to_s}/#{job.end_time.to_s}" : "" , job.user, project_number, job.next_dataset_id]}

    @total = @job_list.length

    # pager
    current_page = params[:format]
    @current_page = (current_page||1).to_i
    @page_list = (1..(@job_list.length.to_f/@page_unit).ceil).to_a
    start = (@current_page - 1) * @page_unit
    last  = @current_page * @page_unit - 1
    @job_list = @job_list[start..last]
  end
  def print_log
    text = 'no log found'
    if @job_id = params[:job_id] and job = Job.find_by_id(@job_id) and
      stdout_path = job.stdout_path and File.exist?(stdout_path) and
      stderr_path = job.stderr_path and File.exist?(stderr_path)
      stdout_text = File.read(stdout_path)
      stderr_text = File.read(stderr_path)
      text = [stdout_path, "-"*50, stdout_text, "___STDOUT_END___\n", stderr_path, "-"*50, stderr_text, "___STDERR_END___"].join("\n")
    end
    render :plain => text
  end
  def print_script
    text = 'no script found'
    if @job_id = params[:job_id] and job = Job.find_by_id(@job_id) and
      script_path = job.script_path and File.exist?(script_path)
      text = File.read(script_path)
    else
      text = @@workflow_manager.get_script(params[:job_id])
    end
    render :plain => text
  end
  def kill_job
    @status = 'kill job failed'
    if @job_id = params[:id] and job = Job.find_by_id(@job_id)
      job.status = "KILLME"
      job.save
      @status = "Killing the job (job ID: #{@job_id})"
      @command = "scancel #{job.submit_job_id}"
    end
  end
  def multi_kill_job
    @job_ids = if flag=params[:kill_flag]
                      flag.keys
                    end
    @statuses = ''
    @commands = ''
    @job_ids.each do |job_id|
      if job = Job.find_by_id(job_id)
        job.status = "KILLME"
        job.save
      end
      @statuses << "Killing the job (job ID: #{job_id})" + "\n"
      @commands << "scancel #{job.submit_job_id}\n"
    end
  end
  def resubmit_job
    if @job_id = params[:id] and job = Job.find_by_id(@job_id)
      @data_set_id = if data_set = job.data_set
                       data_set.id
                     end

      file_base_name = File.basename(job.script_path)
      submit_job_script_dir = SushiFabric::Application.config.submit_job_script_dir
      new_stdout_path = File.join(submit_job_script_dir, file_base_name + "_o.log")
      new_stderr_path = File.join(submit_job_script_dir, file_base_name + "_e.log")
      job.stdout_path = new_stdout_path
      job.stderr_path = new_stderr_path
      job.status = "CREATED"
      job.save

      puts "RESUBMITTED"
    else
      #raise "SOMETHING WRONG"
      puts "FAILED resubmission"
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
