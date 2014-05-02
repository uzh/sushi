class HomeController < ApplicationController
  def index
    @fgcz = if `hostname`.chomp =~ /fgcz-s-034/
              true
            else
              false
            end
    session[:employee] = true if @fgcz and FGCZ.get_user_groups(current_user.login).include?('Employees')
    session[:projects] = if @fgcz 
                           FGCZ.get_user_projects(current_user.login).map{|project| project.gsub(/p/,'').to_i}.sort
                         else
                           [1001]
                         end
    session[:project] = if @fgcz
                          if project=params[:select_project] and number=project[:number] and number.to_i!=0 or
                             project=params[:project] and number=project[:number] and number.to_i!=0 and 
                             (session[:employee] or session[:projects].include?(number.to_i))
                            current_user.selected_project = number
                            current_user.save
                            number.to_i
                          elsif current_user.selected_project != -1
                            current_user.selected_project
                          else
                            session[:projects].first
                          end
                        else
                          session[:projects].first
                        end
    if @fgcz and current_user.selected_project == -1
      current_user.selected_project = session[:project]
      current_user.save
    end
  end
  def gstore
    # path
    @base = '/projects'
    @path = params[:project_id]
    @path = File.join(@path, params[:dirs]) if params[:dirs]
    parent = @path.split('/')
    parent.pop
    @parent = parent.join('/')
    @parent = params[:project_id] if @parent.empty?
    @sort = params[:format] if params[:format]
    @files = Dir[File.join(SushiFabric::GSTORE_DIR, @path)+"/*"]
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
        @files.sort_by! {|file| File.mtime(file)}
      when 'Size'
        @files.sort_by! {|file| File.size(file)}
      end
      @files.reverse! if session[:gstore_reverse]
    else
      session[:gstore_reverse] = nil
      if @path == params[:project_id] 
        @files.sort_by! {|file| File.mtime(file)}
        @files.reverse!
      end
    end
  end
  def sushi_rank
    command = "wfm_job_list -d #{SushiFabric::WORKFLOW_MANAGER}"
    count_name = {}
    first_date = []
    this_month = Time.now.to_s.split.first.split(/-/)[0,2].join('-')
    monthly_mvp = {}
    @count_month = {}
    IO.popen(command) do |io|
      while line=io.gets
        # e.g. "564201,fail,QC_sample_dataset.sh,2013-07-12 17:38:21,masaomi,1001\n"
        id, stat, script, date, name, project = line.chomp.split(/,/)
        count_name[name]||=0
        count_name[name]+=1
        date = date.split.first.split(/-/)[0,2].join('-')
        if date == this_month
          monthly_mvp[name]||=0
          monthly_mvp[name]+=1
        end
        @count_month[date]||=0
        @count_month[date]+=1
        first_date << date
      end
    end
    @count_month = @count_month.to_a.sort
    @rank = count_name.sort_by{|name, count| count}.reverse
    @monthly_mvp = monthly_mvp
    @mvp = monthly_mvp.sort_by{|name, count| count}.reverse.first.first
    @first_date = first_date.sort.first.split.first
  end
end
