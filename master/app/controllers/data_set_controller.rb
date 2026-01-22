class DataSetController < ApplicationController

  skip_before_action :authenticate_user!, only: [:update_completed_samples]

  include SushiFabric
  def top(n_dataset=1000, tree=nil)
    view_context.project_init
    @project = Project.find_by_number(session[:project].to_i)

    @tree = tree
    @data_sets = []
    if @project and  data_sets = @project.data_sets
       @data_sets = data_sets.reverse[0, n_dataset]
    end
  end
  def update_completed_samples_(id)
    sample_available = 0
    if data_set = DataSet.find_by_id(id.to_i)
      sample_available = data_set.samples_length
      if data_set.completed_samples.to_i != data_set.samples_length
        sample_available = 0
        data_set.samples.each do |sample|
          if sample_file = sample.to_hash.select{|header, file| header and header.tag?('File')}.first
            # Check if file value is not nil before splitting
            if sample_file.last && !sample_file.last.to_s.empty?
              file_list = sample_file.last.split(",") ## sample_file is an array holding the header and the file
              all_files_exist = file_list.all? { |f| File.exist?(File.join(SushiFabric::GSTORE_DIR, f)) }
              if all_files_exist
                sample_available+=1
              end
            end
          else # in case of no [File] tag sample
            sample_available+=1
          end
        end
      end
      data_set.completed_samples = sample_available
      data_set.save
    end
    sample_available
  end
  def update_completed_samples
    id = params[:id]
    sample_available = update_completed_samples_(id)
    render plain: sample_available
  end
  def index
    @fgcz = SushiFabric::Application.config.fgcz?
    if warning = session['import_fail']
      @warning = warning
      session['import_fail'] = nil
    end
    top(10)
  end
  def index_full
    top
    render action: "index"
  end
  def index_tree
    top(10, :tree)
    render action: "index"
  end
  def list
    @project = Project.find_by_number(session[:project].to_i)
    @sushi_app = params.dig(:sushi, :app)
    @retired_sushi_app = params.dig(:retired_sushi, :app)
    if @project
      @data_sets = if sushi_app_name = params.dig(:select, :sushi_app) || params.dig(:select, :retired_sushi_app)
                     if sushi_app_name == 'ImportedDataSets'
                      data_sets_ = DataSet.all.select{|data_set|
                        data_set.sushi_app_name.nil?
                      }
                     elsif sushi_app_name == 'AllDataSets'
                      data_sets_ = DataSet.all
                     else
                      data_sets_ = DataSet.all.select{|data_set|
                        data_set.sushi_app_name =~ /#{sushi_app_name}/i
                      }
                     end
                   else
                     data_sets_ = DataSet.all.sort_by{|data_set| Time.now-data_set.created_at}[0,10]
                   end
    end
    @sushi_apps = ["--- select ---", "ImportedDataSets", "AllDataSets"].concat(SushiApplication.all.sort_by{|app| app.class_name}.map{|app| app.class_name})
    lib_dir = File.expand_path('../../../lib/otherApps', __FILE__)
    @retired_sushi_apps = ["--- select ---"].concat(Dir[File.join(lib_dir, '*.rb')].to_a.map{|script| File.basename(script).gsub(/.rb$/, '')}.sort)
    @total = DataSet.select(:id).size
  end
#  caches_action :report
#  caches_page :report
  def bfabric
    @project = Project.find_by_number(session[:project].to_i)
    op = params[:parameters][:bfabric_option]
    pid = Process.fork do
      Process.fork do
        @project.register_bfabric(op)
      end # grand-child process
    end # child process
    Process.waitpid pid

    index
    render action: "index"
  end
  def report
    @project = Project.find_by_number(session[:project].to_i)
    @tree = []
    node_list = {}
    @root = []
    top_nodes = []
    project_dataset_ids = Hash[*(@project.data_sets.map{|data_set| [data_set.id, true]}.flatten)]
    @project.data_sets.each do |data_set|
      report_link = ""
      if i = data_set.headers.index{|header| header.tag?("Link")} 
        report_base = data_set.samples.first.to_hash[data_set.headers[i]]
        base = File.basename(report_base)
        report_link = if data_set.completed_samples.to_i == data_set.samples_length
                        report_url = File.join('/projects', report_base)
                        "<a href='#{report_url}'>#{base}</a>"
                      else 
                        base
                      end
      end
      node = {"id" => data_set.id, 
              "text" => " <a href='/data_set/#{data_set.id}'>"+data_set.name+'</a> '+ data_set.comment.to_s + report_link,
              "children" => []}
      node_list[data_set.id] = node
      if parent = data_set.data_set and project_dataset_ids[parent.id]
        node_list[parent.id]['children'] << node
      else
        top_nodes << node
      end
      if data_set.id == params[:format].to_i
        @root << node
      end
    end
    @root = top_nodes.reverse if @root.empty?
    @tree.concat @root
  end
  def script_log
    data_set = if id = params[:format]
                 DataSet.find_by_id(id)
               end
    @job_list = if data_set
                  data_set.jobs.map{|job| [job.id, File.basename(job.script_path), File.basename(job.stdout_path)]}
                else
                  {}
                end
  end
  def set_runnable_apps(refresh = true)
    if @data_set = DataSet.find_by_id(params[:id]) or (@data_set_id and @data_set = DataSet.find_by_id(@data_set_id))
      sushi_apps = runnable_application(@data_set.headers, refresh)
      @sushi_apps_category = sushi_apps.map{|app| app.analysis_category}.uniq.sort
      @sushi_apps = {}
      @nfcore_apps = []
      sushi_apps.sort_by{|app| app.class_name.to_s}.each do |app|
        # Separate nf-core apps from regular apps
        if app.class_name.to_s =~ /^NfCore/
          @nfcore_apps << app.class_name.to_s
        else
          @sushi_apps[app.analysis_category] ||= []
          @sushi_apps[app.analysis_category] << app.class_name.to_s
        end
      end
      @data_set.runnable_apps = @sushi_apps
      @data_set.nfcore_apps = @nfcore_apps
      saved = @data_set.save
      
      # Debug output
      puts "DEBUG set_runnable_apps: @nfcore_apps count: #{@nfcore_apps.size}"
      puts "DEBUG set_runnable_apps: @data_set.nfcore_apps count after save: #{@data_set.nfcore_apps.size}"
      puts "DEBUG set_runnable_apps: save result: #{saved}"
      $stdout.flush
    end
  end
  def show
    @fgcz = SushiFabric::Application.config.fgcz?
    # for order_id link detection: /data_set/pXXXX/oYYYY, or /data_set/oYYYY
    if params[:project] =~ /^o(\d+)/ or params[:id] =~ /^o(\d+)/ and order_id = $1 and dataset = DataSet.find_by_order_id(order_id.to_i)
      project_number = dataset.project.number
      params[:project] = "p#{project_number}"
      params[:order_id] = order_id
      if !employee? and !user_projects.include?(project_number.to_i)
        redirect_to :controller => "home", :action => "index"
      else
        render action: "show_by_order_id"
      end
    else
      # switch project (from job_monitoring)
      if project = params[:project]
        session[:project] = project.to_i
      end
      if !employee? and data_set_id = params[:id] and data_set = DataSet.find_by_id(data_set_id) and project_number = data_set.project.number and !user_projects.include?(project_number.to_i)
        redirect_to :controller => "home", :action => "index"
      else
        # data_set comment
        if data_set = params[:data_set] and comment = data_set[:comment] and id = data_set[:id]
          data_set = DataSet.find_by_id(id)
          data_set.comment = comment
          data_set.save
        end
        # new data_set name
        if new_data_set = params[:data_set] and name = new_data_set[:name] and id = new_data_set[:id]
          data_set = DataSet.find_by_id(id)
          data_set.name = name
          data_set.save
        end

        # search by RunName and OrderID
        @data_set = DataSet.find_by_id(params[:id])
        unless @data_set
          if project_number = params[:project_id]
            project_number = project_number.gsub(/p/, '').to_i
            session[:project] = project_number
          end
          if data_sets = DataSet.where(run_name_order_id: params[:id])
            if data_sets_ = data_sets.to_a.select{|data_set| data_set.project.number == project_number}
              if @data_set = data_sets_.sort_by{|data_set| data_set.created_at}.first
                params[:id] = @data_set.id
              end
            end
          end
        end

        if @data_set
          @factor_columns = @data_set.headers.select{|header| header.tag?('Factor')}

          session[:latest_data_set_id] = @data_set.id
          # check some properties
          #if session[:employee]
          if employee?
            if parent_dataset = @data_set.data_set
              if @data_set.child == false
                @can_delete_data_files = true
              end
              if parent_dataset.bfabric_id and !@data_set.bfabric_id
                @can_register_bfabric = true
              end
            else
              unless @data_set.bfabric_id
                @can_register_bfabric = true
              end
            end
          end
          # check session[:project]
          unless session[:project] == @data_set.project.number
            session[:project] = @data_set.project.number
            current_user.selected_project = session[:project]
            current_user.save
          end

          # check real data
          @file_exist = {}
          @sample_path = []
          @sample_invalid_name = {}
          sample_count = 0
          if @data_set
            @data_set.samples.each do |sample|
              sample_count+=1
              sample.to_hash.each do |header, file|
                if header and (header.tag?('File') or header.tag?('Link') and file !~ /^http/)
                  if file
                    file_list = file.split(",")
                    @sample_path = @sample_path.concat(file_list.map { |f| File.dirname(f)}.uniq)
                    @file_exist[file] = file_list.all? { |f| File.exist?(File.join(SushiFabric::GSTORE_DIR, f)) }
                  else
                    @file_exist[header] = false
                  end
                else
                  @file_exist[file] = true
                end
                if header == 'Name' and file =~ /[!@\#$%^&*\(\)\<\>\{\}\[\]\/:; '"=+\|]/
                  @sample_invalid_name[file] = true
                end
              end
            end
          end
          @sample_path.uniq!
          @dataset_path = @sample_path.map{|path| path.split('/')[0,2].join('/')}
          @dataset_path.uniq!

          # update num_samples
          if @data_set.num_samples.to_i != sample_count
            @data_set.num_samples = sample_count
          end
          if @data_set.num_samples.to_i != @data_set.completed_samples.to_i
            update_completed_samples_(@data_set.id)
          end

          if !@data_set.refreshed_apps and @data_set.runnable_apps.empty?
            @data_set.refreshed_apps = true
            @data_set.save
            set_runnable_apps(false)
          end

          @employee_apps = employee_apps
          @sushi_apps = @data_set.runnable_apps
          @sushi_apps_category = @sushi_apps.keys.sort
          @nfcore_apps = @data_set.nfcore_apps || []
          
          # Debug output for nf-core apps
          puts "DEBUG show: @sushi_apps keys: #{@sushi_apps.keys.inspect}"
          puts "DEBUG show: @nfcore_apps count: #{@nfcore_apps.size}"
          puts "DEBUG show: @nfcore_apps first 5: #{@nfcore_apps.first(5).inspect}"
          $stdout.flush
        else
          @url_not_found = true
          index
          render action: "index"
        end
      end
    end
  end
  def show_by_order_id
    @fgcz = SushiFabric::Application.config.fgcz?

    order_id = if params[:id] =~ /^o(\d+)/
                 $1
               else
                 params[:id]
               end
    if dataset = DataSet.find_by_order_id(order_id.to_i)
      project_number = dataset.project.number
      params[:project] = "p#{project_number}"
      # params[:id] = nil
      # view_context.project_init
      params[:order_id] = order_id
      #if !session[:employee] and !session[:projects].include?(project_number.to_i)
      if !employee? and !user_projects.include?(project_number.to_i)
        redirect_to :controller => "home", :action => "index"
      end
    end
  end
  def add_comment
    if id = params[:data_set_id] and comment = params[:data_set_comment]
      data_set = DataSet.find_by_id(id)
      data_set.comment = comment
      data_set.save
    end 
    redirect_to(:action => "show") and return
  end
  def edit_name
    if id = params[:data_set_id] and name = params[:data_set_name]
      data_set = DataSet.find_by_id(id)
      data_set.name = name
      data_set.save
      session[:latest_data_set_id] = data_set.id
    end
    redirect_to(:action => "show") and return
  end
  def mod_bfab_data_set_id
    if id = params[:data_set_id] and bfabric_id = params[:bfab_data_set_id]
      data_set = DataSet.find_by_id(id)
      data_set.bfabric_id = bfabric_id
      data_set.save
    end
    redirect_to(:action => "show") and return
  end
  def refresh_apps
    set_runnable_apps
    show
  end
  def edit
    show
  end
  def trace_treeviews(root, data_set, parent_id, project_number, current_data_set, state_opened, data_set_ids={})
    data_set_id = data_set.id
    node_text = if data_set == current_data_set
             "<b>" + data_set.data_sets.length.to_s+" "+data_set.name+"</b> "+" <small><font color='gray'>"+data_set.comment.to_s+"</font></small>"
           else
              data_set.data_sets.length.to_s+" "+data_set.name+" "+" <small><font color='gray'>"+data_set.comment.to_s+"</font></small>"
           end
    node = {"id" => data_set_id, 
            "text" => node_text,
            "parent" => parent_id,
            "state" => {"opened":state_opened},
            "a_attr" => {"href"=>"/data_set/p#{project_number}/#{data_set_id}", 
                         "onclick"=>"location.href = '/data_set/p#{project_number}/#{data_set_id}'"}
            }
    root << node
    data_set_ids[data_set_id] = true
    data_set.data_sets.each do |child|
      if child.project.number==project_number
        trace_treeviews(root, child, data_set.id, project_number, current_data_set, state_opened, data_set_ids)
      end
    end
  end
  def back_trace_treeviews(tree, data_set, data_set_ids={})
    parent_id = if parent = data_set.data_set
                  parent.id
                else
                  "#"
                end
    node_text = data_set.data_sets.length.to_s+" "+data_set.name+" "+" <small><font color='gray'>"+data_set.comment.to_s+"</font></small>"
    data_set_id = data_set.id
    project_number = data_set.project.number
    node = {"id" => data_set_id, 
            "text" => node_text,
            "parent" => parent_id,
            "state" => {"opened":true},
            "a_attr" => {"href"=>"/data_set/p#{project_number}/#{data_set_id}", 
                         "onclick"=>"location.href = '/data_set/p#{project_number}/#{data_set_id}'"}
            }
    tree << node
    data_set_ids[data_set_id] = true
    if parent
      back_trace_treeviews(tree, parent, data_set_ids)
    end
  end
  def partial_treeviews
    root = []
    if current_data_set_id = params[:format]
      # search root parental dataset
      data_set = DataSet.find_by_id(current_data_set_id.to_i)
      parent_id = if parent = data_set.data_set
                    back_trace_treeviews(root, parent)
                    parent.id
                  else
                    "#"
                  end
      state_opened = false
      trace_treeviews(root, data_set, parent_id, data_set.project.number, data_set, state_opened)
    end
    render :json => root.sort_by{|node| node["id"]}.reverse
  end
  def partial_treeviews2
    root = []
    data_set_ids = {}
    if current_data_set_id = params[:format]
      # search root parental dataset
      data_set = DataSet.find_by_id(current_data_set_id.to_i)
      parent_id = if parent = data_set.data_set
                    back_trace_treeviews(root, parent, data_set_ids)
                    parent.id
                  else
                    "#"
                  end
      state_opened = false
      trace_treeviews(root, data_set, parent_id, data_set.project.number, data_set, state_opened, data_set_ids)
    end

    @project = Project.find_by_number(session[:project].to_i)
    project_dataset_ids = Hash[*(@project.data_sets.map{|data_set| [data_set.id, true]}.flatten)]
    @project.data_sets.each do |data_set|
      unless data_set_ids[data_set.id]
        node = {"id" => data_set.id, 
                "text" => data_set.data_sets.length.to_s+" "+data_set.name+" <small><font color='gray'>"+data_set.comment.to_s+"</font></small>",
                "a_attr" => {"href"=>"/data_set/p#{@project.number}/#{data_set.id}", 
                             "onclick"=>"location.href = '/data_set/p#{@project.number}/#{data_set.id}'"}
                }
        if parent = data_set.data_set and project_dataset_ids[parent.id]
          node["parent"] = parent.id
        else
          node["parent"] = "#"
        end
        root << node
      end
    end
 
    render :json => root.sort_by{|node| node["id"]}.reverse
  end
  def partial_treeviews_by_order_id
    tree = []
    if order_id = params[:format] and data_sets = DataSet.select{|ds| ds.order_id == order_id.to_i} and !data_sets.empty?
      data_sets.each do |data_set|
        root = []
        parent_id = "#"
        state_opened = false
        trace_treeviews(root, data_set, parent_id, data_set.project.number, data_set, state_opened)
        tree += root
      end

    end

    render :json => tree.sort_by{|node| node["id"]}.reverse
  end
  def make_whole_tree
    unless @project
      @project = Project.find_by_number(session[:project].to_i)
    end

    project_dataset_ids = @project.data_sets.pluck(:id).index_with(true)
    node_infos = @project.data_sets.includes(:data_sets).map do |data_set|
      [data_set.id, data_set.name, data_set.comment, data_set.parent_id, data_set.data_sets.length]
    end

    project_number = @project.number
    root = node_infos.map do |node_info|
      data_set_id, name, comment, parent_id, num_children = node_info
      {
        "id" => data_set_id,
        "text" => "#{num_children} #{name} <small><font color='gray'>#{comment}</font></small>",
        "a_attr" => { "href" => "/data_set/p#{project_number}/#{data_set_id}" },
        "parent" => (parent_id && project_dataset_ids[parent_id] ? parent_id : "#")
      }
    end
    json = root.sort_by{|node| node["id"]}.reverse.to_json
  end
  def whole_treeviews
    render :json => make_whole_tree
  end
  def import_from_gstore
    params[:project] = session[:project]

    if session[:project] 
      unless @project = Project.find_by_number(session[:project].to_i)
        @project = Project.new
        @project.number = session[:project].to_i
        @project.save
      end

      tsv = File.join(SushiFabric::GSTORE_DIR, "#{params[:dataset]}.#{params[:format]}")
      multi_data_sets = false
      open(tsv) do |input|
        while line=input.gets
          if line =~ /ProjectNumber/
            multi_data_sets = true
            break
          end
        end
      end
      
      if multi_data_sets
        csv = CSV.readlines(tsv, :col_sep=>"\t")
        data_set = []
        headers = []
        rows = []
        csv.each do |row|
          if data_set.empty?
            data_set = row
          elsif headers.empty?
            headers = row
          elsif !row.empty? and !row.join.strip.empty?
            rows << row
          else
            unless headers.include?(nil)
              @data_set_id = DataSet.save_dataset_to_database(
                  data_set_arr: data_set,
                  headers: headers,
                  rows: rows,
                  user: current_user
                )
            else
              session['import_fail'] = 'There must be a blank column. Please check it. Import is incomplete.'
            end
          end
          if row.empty?
            data_set = []
            headers = []
            rows = []
          end
        end
      else
        data_set_tsv = CSV.readlines(tsv, :headers => true, :col_sep=>"\t")

        data_set = []
        headers = data_set_tsv.headers
        rows = []
        items = params[:dataset].split(/\//)
        data_set << "DataSetName"
        data_set << items[-2]
        data_set << "ProjectNumber" << @project.number
        order_ids = {}
        data_set_tsv.each do |row|
          unless row.fields.join.strip.empty?
            rows << row.fields
            if order_id = row["Order Id [B-Fabric]"]
              order_ids[order_id] = true
            end
          end
        end
        unless headers.include?(nil)
          @data_set_id = DataSet.save_dataset_to_database(
              data_set_arr: data_set,
              headers: headers,
              rows: rows,
              user: current_user
            )
        else
          session['import_fail'] = 'There must be a blank column. Please check it. Import is incomplete.'
        end
      end

      if @data_set_id
        refresh = if SushiApplication.count == 0
                    true
                  else
                    false
                  end
        set_runnable_apps(refresh)

        data_set = DataSet.find_by_id(@data_set_id)
        unless order_ids.keys.empty?
          data_set.order_ids = order_ids.keys
          data_set.save
        end

        unless session[:off_bfabric_registration]
          pid = Process.fork do
            Process.fork do
              data_set.register_bfabric
            end # grand-child process
          end # child process
          Process.waitpid pid
        end
        update_completed_samples_(@data_set_id)
      end
    end

    redirect_to :controller => "data_set"
  end
  class CSV::Table
    # imitate DataSet object for call save_dataset_tsv_in_gstore
    def tsv_string
      string = CSV.generate(:col_sep=>"\t") do |out|
        out << headers
        each do |sample|
          out << headers.map{|header|
            val = sample[header]
            val.to_s.empty? ? nil:val}
        end
      end
      string
    end
    def save_as_tsv(file_name)
      File.write(file_name, tsv_string)
    end
    def child
      nil
    end
    def paths
      dirs = []
      each do |sample|
        sample.each do |header, file|
          if (header.tag?('File') or header.tag?('Link')) and file !~ /http/
            dirs << File.dirname(file)
          end
        end
      end
      dirs = dirs.map{|path| path.split('/')[0,2].join('/')}.uniq
    end
    def name
      "DataSet_#{Time.now.strftime("%Y-%m-%d--%H-%M-%S")}"
    end
  end

  def import
    params[:project] = session[:project]

    if session[:project] 
      unless @project = Project.find_by_number(session[:project].to_i)
        @project = Project.new
        @project.number = session[:project].to_i
        @project.save
      end
      @data_set_ids = @project.data_sets.map{|data_set| data_set.id}.push('').reverse

      if file = params[:file] and tsv = file[:name]
        # check meta-dataset
        make_meta_dataset = (meta_dataset = params[:meta_dataset] and plate_samples = meta_dataset[:plate_samples].to_i and plate_samples > 0)
        if make_meta_dataset
          # first, make plate-dataset.tsv
          data_set_tsv = CSV.readlines(tsv.path, :headers => true, :col_sep=>"\t")
          data_set = []
          headers = ["Name", "Species", "CellDataset [File]"]
          rows = []
          data_set << "DataSetName"
          dataset_folder_name = if dataset = params[:dataset] and dataset_name = dataset[:name]
                                  dataset_name
                                else
                                  "MetaDataSet_#{Time.now.strftime("%Y%m%d-%H%M%S")}"
                                end
          data_set << dataset_folder_name
          data_set << "ProjectNumber" << @project.number
          if parent = params[:parent] and parent_id = parent[:id] and parent_data_set = DataSet.find_by_id(parent_id.to_i)
            data_set << "ParentID" << parent_data_set.id
          end
          data_set_first = data_set_tsv.first
          species = (data_set_first["Species"]||"Unknown")
          dataset_path = if tag_file = data_set_first.select{|k,v| k.tag?("File")} and !tag_file.empty?
                           File.dirname(tag_file.first.last)
                         else
                           File.join("p#{@project.number}", dataset_name)
                         end
          plates = data_set_tsv.each_slice(plate_samples).to_a
          digit = ((Math.log10(plates.length)*10).to_i/10)+1
          data_set_tsv.each.with_index do |row, i|
            if dataset_path.nil? and tag_file = row.select{|k,v| k.tag?("File")} and !tag_file.empty?
              dataset_path = File.dirname(tag_file.first.last)
            end
            plate_number = (i / plate_samples) + 1
            if rows.length < plate_number
              plate_name = "Plate_#{"%0#{digit}d" % plate_number}"
              plate_file = "#{plate_name}-dataset.tsv"
              rows << [plate_name, species, File.join(dataset_path, plate_file)]
            end
          end
          unless headers.include?(nil)
            @data_set_id = DataSet.save_dataset_to_database(
                data_set_arr: data_set,
                headers: headers,
                rows: rows,
                user: current_user
              )
          else
            @warning = 'There must be a blank column. Please check it. Import is incomplete.'
          end
          if @data_set_id
            data_set = DataSet.find_by_id(@data_set_id)
            save_dataset_tsv_in_gstore(data_set, "meta-dataset.tsv")
          end
          # second, mata-dataset.tsv
          plates.each.with_index do |plate, i|
            plate_number = i + 1
            plate_name = "Plate_#{"%0#{digit}d" % plate_number}"
            plate_file = "#{plate_name}-dataset.tsv"
            csv_table_plate = CSV::Table.new(plate)
            save_dataset_tsv_in_gstore(csv_table_plate, plate_file)
          end

        else
          data_set_tsv = CSV.readlines(tsv.path, :headers => true, :col_sep=>"\t")

          data_set = []
          headers = data_set_tsv.headers
          rows = []
          order_ids = {}
          data_set << "DataSetName"
          if dataset = params[:dataset] and dataset_name = dataset[:name]
            data_set << dataset_name
          else
            data_set << "DataSet " + (DataSet.all.length+1).to_s
          end
          data_set << "ProjectNumber" << @project.number
          if parent = params[:parent] and parent_id = parent[:id] and parent_data_set = DataSet.find_by_id(parent_id.to_i)
            data_set << "ParentID" << parent_data_set.id
          end
          data_set_tsv.each do |row|
            unless row.fields.join.strip.empty?
              rows << row.fields
              if order_id = row["Order Id [B-Fabric]"]
                order_ids[order_id] = true
              end
            end
          end
          unless headers.include?(nil)
            @data_set_id = DataSet.save_dataset_to_database(
                data_set_arr: data_set,
                headers: headers,
                rows: rows,
                user: current_user
              )
          else
            @warning = 'There must be a blank column. Please check it. Import is incomplete.'
          end
        end
      end

      # ManGO RunName_oBfabricID save
      if @data_set_id and run = params[:run] and run_id = run[:id]
        data_set = DataSet.find_by_id(@data_set_id)
        data_set.run_name_order_id = run_id
        data_set.save
      end

      if @data_set_id
        set_runnable_apps(false)

        data_set = DataSet.find_by_id(@data_set_id)
        unless order_ids.keys.empty?
          data_set.order_ids = order_ids.keys
          data_set.save
        end

        unless session[:off_bfabric_registration]
          pid = Process.fork do
            Process.fork do
              data_set.register_bfabric
            end # grand-child process
          end # child process
          Process.waitpid pid
        end
        update_completed_samples_(@data_set_id)
      elsif file = params[:file] and tsv = file[:name] and @warning.nil?
        @warning = "There might be the same DataSet that has exactly same samples saved in SUSHI. Please check it."
      end
    end
  end
  def save_as_tsv
    tsv_string = 'Error:DataSet is not found'
    data_set_name = if id = params[:id] and data_set = DataSet.find_by_id(id)
                      tsv_string = data_set.tsv_string
                      data_set.name
                    else
                      'dataset'
                    end
     send_data tsv_string,
     :type => 'text/csv',
     :disposition => "attachment; filename=#{data_set_name}.tsv" 
  end
  def data_sets_tsv_string(data_sets)
    headers = ["ID", "Name", "Project", "SushiApp", "Samples", "Who", "Created", "BFabricID"]
    tsv_string = CSV.generate("", headers: headers, write_headers: true, col_sep:"\t") do |out|
      data_sets.each do |data_set|
        row = []
        row << data_set.id
        row << data_set.name
        row << data_set.project.number
        row << data_set.sushi_app_name
        row << "#{data_set.completed_samples.to_i} / #{data_set.samples_length}"
        row << if user = data_set.user
          user.login
        else
          "sushi_lover"
        end
        row << data_set.created_at.strftime("%Y-%b-%d %X ") + SushiFabric::Application.config.time_zone.split('/').last
        row << data_set.bfabric_id.to_s
        out << row
      end
    end
    tsv_string
  end
  def save_all_dataset_list_as_tsv
    project = Project.find_by_number(session[:project].to_i)
    data_sets = if sushi_app_name = params[:sushi_app]
                  if sushi_app_name == 'ImportedDataSets'
                   data_sets_ = DataSet.all.select{|data_set|
                     data_set.sushi_app_name.nil?
                   }
                  elsif sushi_app_name == 'AllDataSets'
                    data_sets_ = DataSet.all
                  else
                   data_sets_ = DataSet.all.select{|data_set|
                     data_set.sushi_app_name =~ /#{sushi_app_name}/i
                   }
                  end
                else
                  data_sets_ = DataSet.all.sort_by{|data_set| Time.now-data_set.created_at}[0,10]
                end
    tsv_string = data_sets_tsv_string(data_sets)
    data_set_name = if sushi_app_name
                      "#{sushi_app_name}_datasets"
                    else
                      "datasets"
                    end
    send_data tsv_string,
    :type => 'text/csv',
    :disposition => "attachment; filename=#{data_set_name}.tsv"
  end
  def save_project_dataset_list_as_tsv
    project = Project.find_by_number(session[:project].to_i)
    data_sets = project.data_sets
    tsv_string = data_sets_tsv_string(data_sets)
    data_set_name = "p#{project.number}_datasets"
    send_data tsv_string,
    :type => 'text/csv',
    :disposition => "attachment; filename=#{data_set_name}.tsv"
  end
  def delete
    @data_set = DataSet.find_by_id(params[:format])

    # check real data
    @file_exist = {}
    @sample_path = []
    @data_set.samples.each do |sample|
      sample.to_hash.each do |header, file|
        if header and file and header.tag?('File')
          file_path = File.join(SushiFabric::GSTORE_DIR, file)
          @sample_path << File.dirname(file)
          @file_exist[file] = File.exist?(file_path)
        else
          @file_exist[file] = true
        end
      end
    end
    @sample_path.uniq!

  end
  def multi_delete
    @data_set_ids = if flag=params[:delete_flag]
                      flag.keys
                    end
    @gstore_dataset_deletable = false
    if @data_set_ids.length == 1
      # same as delete action
      @data_set = DataSet.find_by_id(@data_set_ids.first)

      # check real data
      @file_exist = {}
      @sample_path = []
      @data_set.samples.each do |sample|
        sample.to_hash.each do |header, file|
          if header and file and header.tag?('File')
            file_path = File.join(SushiFabric::GSTORE_DIR, file)
            @sample_path << File.dirname(file)
            @file_exist[file] = File.exist?(file_path)
          else
            @file_exist[file] = true
          end
        end
      end
      @sample_path.uniq!
      if @data_set.parent_id and @data_set.child == false
        @gstore_dataset_deletable = true
      end
      render action: "delete"
    else
      # @data_set_ids.length should_be > 1
      @data_sets = []
      @orig_datasets = []
      @child_datasets = []
      @data_set_ids.each do |id|
        data_set = DataSet.find_by_id(id)
        @data_sets << data_set
        unless data_set.parent_id
          @orig_datasets << data_set
        end
        if data_set.child
          @child_datasets << data_set
        end
      end
      if @orig_datasets.empty? and @child_datasets.empty?
        @gstore_dataset_deletable = true
      end
    end
  end
  def destroy
    @fgcz = SushiFabric::Application.config.fgcz?
    if @data_set = DataSet.find_by_id(params[:id])
      @parent_dataset = @data_set.data_set
      @project_id = @data_set.project.id
      @option = params[:option_delete]

      # check real data
      @sample_path = []
      @data_set.samples.each do |sample|
        sample.to_hash.each do |header, file|
          if header and file and header.tag?('File') 
            file_path = File.join(SushiFabric::GSTORE_DIR, file)
            @sample_path << File.dirname(file)
          end
        end
      end
      @sample_path.uniq!

      # delete data in gstore
      if @sample_path.first
        target = File.join(SushiFabric::GSTORE_DIR, @sample_path.first)
        @command = @@sushi_server.delete_command(target, bfabric_dataset_id: @data_set.bfabric_id)
        if @option[:only_gstore] == "1"
          @command_log = `#{@command}`
          if request = @command_log.split and request_no = request[4]
            @greq_status_command = "g-req status #{request_no}"
          end
        end
      end

      # delete data in sushi
      if @option[:only_sushi] == "1"
        @data_set.samples.each do |sample|
          sample.delete
        end
        @deleted_data_set = @data_set.delete
      else
        @deleted_data_set = @data_set
      end

      # delete data files
      if @option[:data_files] == '1'
        @data_set = DataSet.find_by_id(params[:id])
        delete_candidates(@data_set)
        render action: "confirm_delete_only_data_files"
      end
      params[:project_id] = "p#{session[:project]}" if session[:project]
      @deleted_data_set
    end
  end
  def multi_destroy
    @option = params[:option_delete]
    @data_set_ids = params[:option][:data_set_ids].split(',')
    @commands = []
    @command_logs = []
    @deleted_data_sets = []
    @parent_datasets = @data_set_ids.map{|id|
      if data_set = DataSet.find_by_id(id)
        parent_dataset = data_set.data_set
      end
    }.compact.uniq.sort_by{|dataset| dataset.id}

    @data_set_ids.each do |id|
      params[:id] = id
      @deleted_data_sets << destroy
      if @command
        @commands << @command.chomp
        @command = nil
      end
      if @command_log
        @command_logs << @command_log.chomp
        @command_log = nil
      end
    end
    params[:project_id] = "p#{session[:project]}" if session[:project]
  end
  def job_parameter
    @data_set = if id = params[:format]
                  DataSet.find_by_id(id)
                end
    if @parameters = @data_set.job_parameters and @parameters.empty?
      @sample_path = if @data_set
                       sample_path(@data_set)
                     end
      @parameters_tsv = if @sample_path
                          File.join(SushiFabric::GSTORE_DIR, @sample_path, 'parameters.tsv')
                        end
      if @parameters_tsv and File.exist?(@parameters_tsv)
        File.readlines(@parameters_tsv).each do |line|
          header, *values = line.chomp.split
          @parameters[header] = values.join(" ")
        end
      end
    end
  end
  def project_paths(data_set)
    paths = []
    if sample = data_set.samples.first
      sample.to_hash.each do |header, file|
        if header and (header.tag?('File') or header.tag?('Link') and file !~ /^http/) and file
          project_path = file.split('/')[0,3].join('/')
          file_path = File.join(SushiFabric::GSTORE_DIR, project_path)
          paths << File.dirname(file_path)
        end
      end
    end
    paths.uniq!
    paths
  end
  def delete_candidates(data_set)
    @delete_candidates = []
    project_paths(data_set).each do |dir|
      Dir[File.join(dir, "*.*")].sort.each do |file|
        unless file =~ /.tsv/ or File.ftype(file) == "directory"
          @delete_candidates << file
        end
      end
    end
    @delete_candidates
  end
  def confirm_delete_only_data_files
    @data_set = DataSet.find_by_id(params[:id])
    delete_candidates(@data_set)
  end
  def run_delete_only_data_files
    @data_set = DataSet.find_by_id(params[:id])
    @delete_files = delete_candidates(@data_set)
    file_exts = @delete_files.map{|file| File.join(File.dirname(file), "*." + file.split('.').last)}.uniq.sort
    target = file_exts.join(" ")
    @command = @@sushi_server.delete_command(target, bfabric_dataset_id: @data_set.bfabric_id)
    @command_log = `#{@command}`
  end
  def register_bfabric
    data_set = DataSet.find_by_id(params[:id])
    # Abort if duplicate headers exist when ignoring tag suffixes (e.g., "[File]", "[Link]")
    if data_set && (dup = data_set.send(:duplicate_headers_ignoring_tags)).any?
      flash[:alert] = "BFabric registration aborted: duplicate column names (ignoring tags): #{dup.join(', ')}"
      redirect_to :controller => "data_set", :action => "show", :id => data_set.id and return
    end
    # Abort if there is no [File]-tagged column
    if data_set && !data_set.send(:has_file_tag_column?)
      flash[:alert] = "BFabric registration aborted: dataset has no [File] tag column"
      redirect_to :controller => "data_set", :action => "show", :id => data_set.id and return
    end
    pid = Process.fork do
      Process.fork do
        #data_set.register_bfabric("only_one")
        data_set.register_bfabric("only_one", register_child_dataset_too: true)
      end # grand-child process
    end # child process
    Process.waitpid pid
    redirect_to :controller => "data_set", :action => "show"
  end
  def update_resource_size
    data_set = DataSet.find_by_id(params[:id])
    pid = Process.fork do
      Process.fork do
        data_set.update_resource_size
      end # grand-child process
    end # child process
    Process.waitpid pid
    redirect_to :controller => "data_set", :action => "show"
  end
  def announce_template_set
    @data_set_id = params[:id]
    @announce_templates = Dir["announce_templates/*.txt"].to_a.sort
  end
  def announce_replace_set
    @template_path = if template = params[:template]
                       template[:path]
                     end
    id = if data_set = params[:data_set]
           data_set[:id]
         end
    @data_set = DataSet.find_by_id(id)
    fastqc_data_set = @data_set.data_sets.select{|dataset| dataset.name =~ /Fastqc/i}.first
    @fastqc_link = if fastqc_data_set and sample = fastqc_data_set.samples.first
                     sample.to_hash["Html [Link]"]
                   else
                     "FASTQC_LINK"
                   end
    fastqscreen_data_set = @data_set.data_sets.select{|dataset| dataset.name =~ /Fastqscreen/i}.first
    @fastqscreen_link = if fastqscreen_data_set and sample = fastqscreen_data_set.samples.first
                     sample.to_hash["Html [Link]"]
                   else
                     "FASTQSCREEN_LINK"
                   end
    @replaces = {}
    @template = []
    File.readlines(@template_path).each do |line|
      @template << line
      if matches = line.scan(/[A-Z_]{2,}/).reject{|x| x=="FGCZ" or x=="SUSHI"}
        matches.each do |key|
          case key
          when "USER_NAME"
            @replaces[key] = "Project members"
          when "PROJECT_NUMBER"
            @replaces[key] = session[:project]
          when "ORDER_NUMBER"
          @replaces[key] = if @data_set.name =~ /_o(\d+)/ or @data_set.name =~ /^o(\d+)_/
                               $1
                             else
                               key
                             end
          when "DATASET_ID"
            @replaces[key] = @data_set.id
          when "DATASET_NAME"
            @replaces[key] = @data_set.name
          when "WORKUNIT_ID"
            @replaces[key] = (@data_set.workunit_id||key)
          when "MY_NAME"
            @replaces[key] = current_user.login.capitalize
          when "FASTQC_LINK"
            @replaces[key] = @fastqc_link.to_s
          when "FASTQSCREEN_LINK"
            @replaces[key] = @fastqscreen_link.to_s
          else
            @replaces[key] = key
          end
        end
      end
    end
  end
  def announce
    params.permit!
    template_path = if template = params[:template]
                      template[:path]
                    end
    @replaces = params[:replaces]
    @template = []
    File.readlines(template_path).each do |line|
      @template << line.chomp.gsub(/#{@replaces.keys.join("|")}/, @replaces)
    end
    @bfab_order_number = if @replaces["DATASET_NAME"] =~ /_o(\d+)/ or @replaces["DATASET_NAME"] =~ /^o(\d+)_/
                           $1
                         end
  end
  
  # Display merge dialog (AJAX)
  def merge_dialog
    @current_dataset = DataSet.find_by_id(params[:id])
    unless @current_dataset
      render plain: "Dataset not found", status: :not_found
      return
    end
    
    @project = @current_dataset.project
    @available_datasets = @project.data_sets.where.not(id: @current_dataset.id).order(created_at: :desc)
    render partial: 'merge_dialog', layout: false
  end
  
  # Merge datasets action
  def merge_with_dataset
    dataset1 = DataSet.find_by_id(params[:id])
    dataset2 = DataSet.find_by_id(params[:second_dataset_id])
    
    unless dataset1 && dataset2
      flash[:error] = "One or both datasets not found"
      redirect_to data_set_path(params[:id]) and return
    end
    
    begin
      merged_dataset = dataset1.merge_with(dataset2, options: {
        merged_dataset_name: params[:merged_name],
        user: current_user,
        excluded_columns: ['Sample Id [B-Fabric]']
      })
      
      flash[:notice] = "Successfully merged datasets. New dataset: #{merged_dataset.name} (ID: #{merged_dataset.id})"
      redirect_to data_set_path(merged_dataset)
    rescue => e
      flash[:error] = "Failed to merge datasets: #{e.message}"
      redirect_to data_set_path(dataset1)
    end
  end
end
