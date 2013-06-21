class RunApplicationController < ApplicationController
  def index
    @sushi_apps = Dir['lib/*.rb'].sort.select{|script| script !~ /sushiApp/ and script !~ /sushiToolBox/ and script !~ /optparse/}.to_a.map{|script| File.basename(script).gsub(/\.rb/,'')}
    if project_number = session[:project] and project = Project.find_by_number(project_number.to_i)
      @data_sets = project.data_sets
    end
  end
  def select_application
    if data_set = params[:data_set] and id = data_set[:id]
      @data_set = DataSet.find_by_id(id.to_i)
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
  def set_parameters
    class_name = params[:app]
    require class_name
    @sushi_app = eval(class_name).new
    data_set_id = params[:data_set][:id]
    @data_set = DataSet.find(data_set_id.to_i)
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
  end
end
