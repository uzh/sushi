class HomeController < ApplicationController
  # Rank reset date - change this value to reset rank aggregation period
  # Past data is preserved in the database; only display/aggregation is affected
  RANK_RESET_DATE = Time.new(2026, 1, 1)
  def index
    @fgcz = SushiFabric::Application.config.fgcz?
    if project_ = params[:project] and project_[:number] =~ /^o(\d+)/ and order_id = $1 and dataset = DataSet.find_by_order_id(order_id.to_i)
      project_number = dataset.project.number
      if session[:employee] or session[:projects].include?(project_number.to_i)
        redirect_to "/data_set/p#{project_number}/o#{order_id}"
      end
    end
    view_context.project_init
  end
  def switch_bfabric_registration
    if bf = params[:bfabric_registration]
      session[:off_bfabric_registration] = (bf == "OFF")
    end
    render action: "index"
  end
  def valid_file_path?(file_path)
    file_path.match?(/\Ap\d{4,}\/[\w\-\._\/]+\z/)
  end
  def gstore
    view_context.project_init
    if !session[:employee] and project_id = params[:project_id] and number = project_id.gsub(/^p/, '') and !session[:projects].include?(number.to_i)
      redirect_to :controller => "home", :action => "index"
    else
      # path
      @base = '/projects'
      @path = params[:project_id]
      @path = File.join(@path, params[:dirs]) if params[:dirs]
      parent = @path.split('/')
      parent.pop
      @parent = parent.join('/')
      @parent = params[:project_id] if @parent.empty?
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
      current_page, @sort, non_rev = params[:format].to_s.split(/:/)
      @current_page = (current_page||1).to_i
      @page_list = (1..(@files.length.to_f/@page_unit).ceil).to_a
      start = (@current_page - 1) * @page_unit
      last  = @current_page * @page_unit - 1

      # sort
      if @sort
        unless non_rev
          session[:gstore_reverse] = !session[:gstore_reverse]
        end
        case @sort
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
      @files = @files[start..last]
      @fgcz = SushiFabric::Application.config.fgcz?

      if @fgcz and !@files
        file_ext = params[:format]
        file_path = "#{@path}.#{file_ext}"
        file_full_path = File.join(SushiFabric::GSTORE_DIR, file_path)

        case file_ext
        when /gz$|bam$|html$/
          if valid_file_path?(file_path)
            redirect_to "https://fgcz-gstore.uzh.ch/projects/#{file_path}", allow_other_host: true
          else
            render plain: "Invalid file path", status: :bad_request
          end
        #when /html$/
        #  send_file file_full_path, disposition: 'inline', type: 'text/html'
        when /log$|tsv$|txt$|sh$/
          send_file file_full_path, disposition: 'inline', type: 'text/plain'
        else
          send_file file_full_path, disposition: 'inline'
        end
      end
    end
  end
  def sushi_rank
    start_of_month = Time.current.beginning_of_month
    end_of_month = Time.current.end_of_month

    # Base scope: only count jobs since the reset date
    rank_jobs = Job.where.not(user: nil).where("created_at >= ?", RANK_RESET_DATE)

    @rank = rank_jobs.group(:user).count.sort_by{|name, count| -count}
    @count_month = rank_jobs
                     .group("DATE_FORMAT(created_at, '%Y-%m')")
                     .count
                     .sort
    @monthly_mvp = rank_jobs
                     .where(created_at: start_of_month..end_of_month)
                     .group(:user)
                     .count
    @mvp = unless @monthly_mvp.empty?
             @monthly_mvp.sort_by{|name, count| count}.reverse.first.first
           else
             {}
           end
    @first_date = @count_month.empty? ? RANK_RESET_DATE.strftime('%Y-%m') : @count_month.first.first
    @rank_reset_date = RANK_RESET_DATE

    subquery = Job.select("user, MIN(created_at) AS first_job_date")
                 .where.not(user: nil)
                 .where("created_at >= ?", RANK_RESET_DATE)
                 .group(:user)
    @count_users = Job.from(subquery, :first_jobs)
                     .group("DATE_FORMAT(first_jobs.first_job_date, '%Y-%m')")
                     .count
                     .sort
    @total_users = []
    total_users = 0
    @count_users.each do |date,users|
      total_users += users.to_i
      @total_users << [date, total_users]
    end
  end
end
