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
    public_dir = File.expand_path('../../../public', __FILE__)
    text = @@workflow_manager.get_log(params[:job_id], :with_err)
    render :plain => text
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
    render :plain => text
  end
  def kill_job
    @status = 'kill job failed'
    if @job_id = params[:id]
      public_dir = File.expand_path('../../../public', __FILE__)
      @status = @@workflow_manager.kill_job(@job_id)
      @command = "wfm_kill_job -i #{@job_id} -d #{SushiFabric::WORKFLOW_MANAGER}"
    end
  end
  def multi_kill_job
    @job_ids = if flag=params[:kill_flag]
                      flag.keys
                    end
    @statuses = ''
    @commands = ''
    @job_ids.each do |job_id|
      @statuses << @@workflow_manager.kill_job(job_id) + "\n"
      @commands << "wfm_kill_job -i #{job_id} -d #{SushiFabric::WORKFLOW_MANAGER}\n"
    end
  end
  def resubmit_job
    if @job_id = params[:id]
      gstore_script_dir = if job = Job.find_by_submit_job_id(@job_id)
                            @data_set_id = if data_set = job.data_set
                                             data_set.id
                                           end
                            File.dirname(job.script_path)
                          end
      prev_params = {}
      if parameters_tsv = File.join(File.dirname(gstore_script_dir), "parameters.tsv") and File.exist?(parameters_tsv)
        File.readlines(parameters_tsv).each do |line|
          name, value = line.chomp.split
          prev_params[name] = value.to_s.delete('"')
        end
        script_content = @@workflow_manager.get_script(@job_id)
        script_path = @@workflow_manager.get_script_path(@job_id)
        project_number = session[:project]
        gsub_options = []
        gsub_options << "-c #{prev_params['cores']}" unless prev_params['cores'].to_s.empty?
        gsub_options << "-n #{prev_params['node']}" unless prev_params['node'].to_s.empty?
        gsub_options << "-p #{prev_params['partition']}" unless prev_params['partition'].to_s.empty?
        gsub_options << "-r #{prev_params['ram']}" unless prev_params['ram'].to_s.empty?
        gsub_options << "-s #{prev_params['scratch']}" unless prev_params['scratch'].to_s.empty?
        gsub_options << "-i #{prev_params['nice']}" unless prev_params['nice'].to_s.empty?
      end

      if script_path and current_user and script_content and project_number and gstore_script_dir and @data_set_id and File.exist?(parameters_tsv)
        new_job_id = @@workflow_manager.start_monitoring3(script_path, script_content, current_user.login, project_number, gsub_options.join(' '), gstore_script_dir, job.next_dataset_id, SushiFabric::RAILS_HOST)
        puts "job_id: #{new_job_id}"
        new_job = Job.new
        new_job.submit_job_id = new_job_id.to_i
        new_job.script_path = script_path
        new_job.next_dataset_id = job.next_dataset_id
        new_job.save
        new_job.data_set.jobs << new_job
        new_job.data_set.save

        puts "RESUBMITTED"
      else
        #raise "SOMETHING WRONG"
        puts "FAILED resubmission"
      end
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
