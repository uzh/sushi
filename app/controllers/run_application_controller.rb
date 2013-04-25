class RunApplicationController < ApplicationController
  def index
    @sushi_apps = Dir['lib/*.rb'].sort.select{|script| script !~ /sushiApp\.rb/}.to_a.map{|script| File.basename(script).gsub(/\.rb/,'')}
    if project_number = session[:project] and project = Project.find_by_number(project_number.to_i)
      @data_sets = project.data_sets
    end
  end
  def set_parameters
    @params = params
    class_name = params[:sushi_app][:class]
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
      @sushi_app.params[key].value = if @sushi_app.params[key].type != String
                                       eval(value)
                                     else
                                       value
                                     end
    end
  end
end
