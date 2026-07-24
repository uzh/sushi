class JobMonitoringController < ApplicationController
  JOB_LIST_COLUMNS =
    'jobs.id, jobs.status, jobs.script_path, jobs.start_time, jobs.end_time, ' \
    'jobs.user, jobs.next_dataset_id, projects.number AS project_number, ' \
    'data_sets.name AS next_dataset_name'.freeze

  def fetch_jobs(params, project_number)
    option = params[:option]

    all_projects =
      if option && option[:all_job_list]        # all projects (button)
        true
      elsif option && option[:project_job_list] # project specific (button)
        false
      else                                       # fall back to session
        !!session[:all_job_list]
      end

    if all_projects
      configure_job_list(all_job_list: true, page_unit: 1000)
      # The all-projects list used to time out (504). With a plain join the
      # optimizer drove from `projects` (a tiny table), so ORDER BY jobs.id could
      # not use the primary key and the entire jobs->data_sets->projects join was
      # materialized and filesorted before LIMIT. STRAIGHT_JOIN forces `jobs` as
      # the driving table, so LIMIT 1000 is served by a backward primary-key scan
      # that stops early. The result set is identical to the old query (the newest
      # 1000 jobs that resolve to a data_set + project).
      Job.select(JOB_LIST_COLUMNS)
         .joins('STRAIGHT_JOIN data_sets ON data_sets.id = jobs.next_dataset_id')
         .joins('STRAIGHT_JOIN projects ON projects.id = data_sets.project_id')
         .order('jobs.id DESC')
         .limit(1000)
    else
      # Project-specific list filters on projects.number first (a small subset),
      # so the optimizer's own join order is fine here; keep the plain join.
      configure_job_list(all_job_list: false)
      Job.select(JOB_LIST_COLUMNS)
         .joins(data_set: :project)
         .where(projects: { number: project_number })
         .order('jobs.id DESC')
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
