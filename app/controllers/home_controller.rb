class HomeController < ApplicationController
  def index
    # tentatively only for develop
    session[:project] = 1001
  end
  def result
    @path = params[:project_id]
    @path = File.join(@path, params[:dirs]) if params[:dirs]
    parent = @path.split('/')
    parent.pop
    @parent = parent.join('/')
    @parent = params[:project_id] if @parent.empty?
  end
end
