class HomeController < ApplicationController
  def index
    # tentatively only for develop
    session[:project] = 1001
  end
  def gstore
    # path
    @base = '/gstore/sushi'
    @path = params[:project_id]
    @path = File.join(@path, params[:dirs]) if params[:dirs]
    parent = @path.split('/')
    parent.pop
    @parent = parent.join('/')
    @parent = params[:project_id] if @parent.empty?
    @sort = params[:format] if params[:format]
    @files = Dir[File.join(GSTORE_DIR, @path)+"/*"]
    @total = @files.length

    # pager
    @page_unit = if page = params[:page] and unit = page[:unit] 
                   session[:gstore_page_unit] = unit.to_i
                 elsif unit = session[:gstore_page_unit]
                   unit.to_i
                 else
                   session[:gstore_page_unit] = 10
                 end
    current_page, sort = params[:format].to_s.split(/:/)
    @current_page = (current_page||1).to_i
    @page_list = (1..(@files.length.to_f/@page_unit).ceil).to_a
    start = (@current_page - 1) * @page_unit
    last  = @current_page * @page_unit - 1
    @files = @files[start..last]

    # sort
    if sort
      session[:gstore_reverse] = !session[:gstore_reverse]
      case sort
      when 'Name'
        @files.sort_by! {|file| File.basename(file)}
      when 'Last_Modified'
        @files.sort_by! {|file| File.ctime(file)}
      when 'Size'
        @files.sort_by! {|file| File.size(file)}
      end
      @files.reverse! if session[:gstore_reverse]
    else
      session[:gstore_reverse] = nil
      if @path == params[:project_id] 
        @files.sort_by! {|file| File.ctime(file)}
        @files.reverse!
      end
    end

  end
end
