class HomeController < ApplicationController
  def index
    # tentatively only for develop
    session[:project] = 1001
  end
  def gstore
    @base = '/gstore/sushi'
    @path = params[:project_id]
    @path = File.join(@path, params[:dirs]) if params[:dirs]
    parent = @path.split('/')
    parent.pop
    @parent = parent.join('/')
    @parent = params[:project_id] if @parent.empty?
    @sort = params[:format] if params[:format]
    @files = Dir[File.join(GSTORE_DIR, @path)+"/*"]
    if params[:format]
      session[:gstore_reverse] = !session[:gstore_reverse]
    else
      session[:gstore_reverse] = nil
    end
    case params[:format]
    when 'Name'
      @files.sort_by! {|file| File.basename(file)}
    when 'Last_Modified'
      @files.sort_by! {|file| File.ctime(file)}
    when 'Size'
      @files.sort_by! {|file| File.size(file)}
    end
    if @path == params[:project_id] and !params[:format]
      @files.sort_by! {|file| File.ctime(file)}
      @files.reverse!
    end
    @files.reverse! if session[:gstore_reverse]
  end
end
