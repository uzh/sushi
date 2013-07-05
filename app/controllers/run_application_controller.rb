class RunApplicationController < ApplicationController
  def index
    @data_sets = if project_number = session[:project] and project = Project.find_by_number(project_number.to_i)
                   project.data_sets
                 else
                   []
                 end
  end
  def select_application
    if data_set_id = params[:format]
      @data_set = DataSet.find_by_id(data_set_id.to_i)

      # check real data
      @file_exist = {}
      @data_set.samples.each do |sample|
        sample.to_hash.values.each do |file|
          if file.split(/\//).first =~ /^p\d+/
            file_path = File.join(GSTORE_DIR, file)
            @file_exist[file] = File.exist?(file_path)
          else
            @file_exist[file] = true
          end
        end
      end

      if @file_exist.values.inject{|a,b| a and b}
        # prepare application buttons
        @sushi_apps = Dir['lib/*.rb'].sort.select{|script| script !~ /sushiApp/ and script !~ /sushiToolBox/ and script !~ /optparse/}.to_a.map{|script| File.basename(script).gsub(/\.rb/,'')}

        # filter application with data_set
        headers = @data_set.headers 
        @sushi_apps = @sushi_apps.select do |class_name|
          require class_name
          sushi_app = eval(class_name).new
          required_columns = sushi_app.required_columns
          (required_columns - headers).empty?
        end
      end
    end
  end
  def set_parameters
    class_name = params[:app]
    require class_name
    @sushi_app = eval(class_name).new
    data_set_id = params[:data_set][:id]
    @data_set = DataSet.find(data_set_id.to_i)
    @nodes = {
      '' => '',
      'fgcz-c-046: cpu 64,mem 504 GB,scr  11T' => 'fgcz-c-046',
      'fgcz-c-047: cpu 32,mem   1 TB,scr  28T' => 'fgcz-c-047',
      'fgcz-c-048: cpu 48,mem 252 GB,scr 3.5T' => 'fgcz-c-048',
      'fgcz-c-049: cpu  8,mem  63 GB,scr 1.7T' => 'fgcz-c-049',
      'fgcz-c-050: cpu  2,mem   3 GB,scr  20G' => 'fgcz-c-050',
      'fgcz-c-051: cpu  8,mem  31 GB,scr 800G' => 'fgcz-c-051',
      'fgcz-c-052: cpu  8,mem  31 GB,scr 800G' => 'fgcz-c-052',
      'fgcz-c-053: cpu  8,mem  31 GB,scr 800G' => 'fgcz-c-053',
      'fgcz-c-054: cpu  8,mem  31 GB,scr 800G' => 'fgcz-c-054',
      'fgcz-c-055: cpu  8,mem  31 GB,scr 800G' => 'fgcz-c-055',
      'fgcz-c-057: cpu  8,mem  31 GB,scr 200G' => 'fgcz-c-057',
      'fgcz-c-058: cpu  8,mem  31 GB,scr 200G' => 'fgcz-c-058',
      'fgcz-c-059: cpu  8,mem  31 GB,scr 200G' => 'fgcz-c-059',
      'fgcz-c-061: cpu  8,mem  31 GB,scr 200G' => 'fgcz-c-061',
      'fgcz-c-063: cpu 12,mem  70 GB,scr 450G' => 'fgcz-c-063',
      'fgcz-c-064: cpu 24,mem  35 GB,scr 4.9T' => 'fgcz-c-064',
      'fgcz-c-065: cpu 24,mem  70 GB,scr 197G' => 'fgcz-c-065',
    }
  end
  def confirmation
    @params = params
    class_name = params[:sushi_app][:class]
    @sushi_app = eval(class_name).new
    data_set_id = params[:data_set][:id]
    @data_set = DataSet.find(data_set_id.to_i)
    params[:parameters].each do |key, value|
      @sushi_app.params[key] = if @sushi_app.params.data_type(key) == String
                                 value
                               else
                                 eval(value)
                               end
    end
  end
  def submit_jobs
    @params = params
    class_name = params[:sushi_app][:class]
    @sushi_app = eval(class_name).new
    @sushi_app.user = current_user.login
    data_set_id = params[:data_set][:id]
    params[:parameters].each do |key, value|
      @sushi_app.params[key] = if @sushi_app.params.data_type(key) == String
                                       value
                                     else
                                       eval(value)
                                     end
    end
    if project_number = session[:project] and project = Project.find_by_number(project_number.to_i)
      @sushi_app.project = 'p' + project_number.to_s
    end
    @sushi_app.dataset_sushi_id = data_set_id.to_i
    @sushi_app.run
    @sushi_app.job_ids.each do |job_id|
      new_job = Job.new
      new_job.submit_job_id = job_id.to_i
      new_job.next_dataset_id = @sushi_app.next_dataset_id
      new_job.save
      new_job.data_set.jobs << new_job
      new_job.data_set.save
    end
  end
end
