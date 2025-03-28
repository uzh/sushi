class JobMonitoringController < ApplicationController
  def fetch_jobs(params, project_number)
    option = params[:option]

    if option && option[:all_job_list] # all projects
      configure_job_list(all_job_list: true, page_unit: 1000)
      Job.joins(data_set: :project)
         .joins('LEFT JOIN data_sets AS next_data_set ON jobs.next_dataset_id = next_data_set.id')
         .select('jobs.*, projects.number AS project_number, next_data_set.name AS next_dataset_name')
         .order(id: :desc)
         .limit(1000)
    elsif option && option[:project_job_list] # project specific
      configure_job_list(all_job_list: false)
      Job.joins(data_set: :project)
         .joins('LEFT JOIN data_sets AS next_data_set ON jobs.next_dataset_id = next_data_set.id')
         .select('jobs.*, projects.number AS project_number, next_data_set.name AS next_dataset_name')
         .where(projects: { number: project_number })
         .order(id: :desc)
         .limit(100)
    elsif session[:all_job_list] # all projects from session
      configure_job_list(all_job_list: true, page_unit: 1000)
      Job.joins(data_set: :project)
         .joins('LEFT JOIN data_sets AS next_data_set ON jobs.next_dataset_id = next_data_set.id')
         .select('jobs.*, projects.number AS project_number, next_data_set.name AS next_dataset_name')
         .order(id: :desc)
         .limit(1000)
    else # project specific
      configure_job_list(all_job_list: false)
      Job.joins(data_set: :project)
         .joins('LEFT JOIN data_sets AS next_data_set ON jobs.next_dataset_id = next_data_set.id')
         .select('jobs.*, projects.number AS project_number, next_data_set.name AS next_dataset_name')
         .where(projects: { number: project_number })
         .order(id: :desc)
         .limit(100)
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

    @job_list = jobs.map{|job|
      start_time = if job.start_time
                     job.start_time.strftime('%Y-%m-%d %H:%M:%S')
                   else
                     ""
                   end
      end_time = if job.end_time
                   job.end_time.strftime('%Y-%m-%d %H:%M:%S')
                 else
                   ""
                 end
      job_script = if job.script_path
                     File.basename(job.script_path)
                   else
                     ""
                   end
      [job.id, job.status.to_s, job_script, job.start_time ? "#{start_time}/#{end_time}" : "" , job.user, job.project_number, job.next_dataset_id, job.next_dataset_name]}

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
      script_text = File.read(script_path)
      text_lines = script_text.lines
      if text_lines.size > 1
        text_lines.insert(1, "##{script_path}")
        text = text_lines.join
      end
    end
    render :plain => text
  end
  def kill_job
    @status = 'kill job failed'
    if @job_id = params[:id] and job = Job.find_by_id(@job_id)
      job.status = "KILL_ME"
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
        job.status = "KILL_ME"
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

      job.status = "CREATED"
      job.save

      puts "RESUBMITTED"
    else
      #raise "SOMETHING WRONG"
      puts "FAILED resubmission"
    end
  end
end
