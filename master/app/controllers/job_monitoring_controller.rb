class JobMonitoringController < ApplicationController
  def index
    public_dir = File.expand_path('../../../public', __FILE__)
    @page_unit = 100
    workflow_manager = DRbObject.new_with_uri(SushiFabric::WORKFLOW_MANAGER)
    @job_list = if option=params[:option] and option[:all_job_list] 
                  @page_unit = 1000
                  @all_job_list=true
                  session[:all_job_list] = true
                  workflow_manager.job_list(false, nil)
                elsif option=params[:option] and option[:project_job_list] 
                  @all_job_list=false
                  session[:all_job_list] = false
                  workflow_manager.job_list(false, session[:project])
                elsif session[:all_job_list]
                  @page_unit = 1000
                  @all_job_list=true
                  session[:all_job_list] = true
                  workflow_manager.job_list(false, nil)
                else
                  @all_job_list=false
                  session[:all_job_list] = false
                  workflow_manager.job_list(false, session[:project])
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
                            data_set = job.data_set
                            @data_set_id = data_set.id
                            File.dirname(job.script_path)
                          end
      parameters_tsv = File.join(File.dirname(gstore_script_dir), "parameters.tsv")
      prev_params = {}
      File.readlines(parameters_tsv).each do |line|
        name, value = line.chomp.split
        prev_params[name] = value.delete('"')
      end 
      script_content = @@workflow_manager.get_script(@job_id)
      script_path = @@workflow_manager.get_script_path(@job_id)
      project_number = session[:project]
      gsub_options = []
      gsub_options << "-c #{prev_params['cores']}" unless prev_params['cores'].to_s.empty?
      gsub_options << "-n #{prev_params['node']}" unless prev_params['node'].to_s.empty?
      gsub_options << "-r #{prev_params['ram']}" unless prev_params['ram'].to_s.empty?
      gsub_options << "-s #{prev_params['scratch']}" unless prev_params['scratch'].to_s.empty?
      if script_path and current_user and script_content and project_number and gstore_script_dir and @data_set_id and File.exist?(parameters_tsv)
        new_job_id = @@workflow_manager.start_monitoring(script_path, current_user.login, 0, script_content, project_number, gsub_options.join(' '), gstore_script_dir)
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
        raise "SOMETHING WRONG"
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
