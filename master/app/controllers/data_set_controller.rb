class DataSetController < ApplicationController
  include SushiFabric
  def top(n_dataset=1000)
    view_context.project_init
    @project = Project.find_by_number(session[:project].to_i)

    @data_sets = []
    if @project and  data_sets = @project.data_sets
      @data_sets = data_sets.reverse[0, n_dataset]
      @data_sets.each do |data_set|
        unless data_set.completed_samples.to_i == data_set.samples_length.to_i
          sample_available = 0
          data_set.samples.each do |sample|
            if sample_file = sample.to_hash.select{|header, file| header and header.tag?('File')}.first
              file_path = File.join(SushiFabric::GSTORE_DIR, sample_file.last) 
              if File.exist?(file_path)
                sample_available+=1
              end
            end
          end
          data_set.completed_samples = sample_available
          data_set.save
        end
      end
    end
  end
  def index
    if warning = session['import_fail']
      @warning = warning
      session['import_fail'] = nil
    end
    top(20)
  end
  def index_full
    top
    render action: "index"
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
    @data_set = if id = params[:format]
                  DataSet.find_by_id(id)
                end
  end
  def set_runnable_apps(refresh = true)
    if @data_set = DataSet.find_by_id(params[:id]) or (@data_set_id and @data_set = DataSet.find_by_id(@data_set_id))
      sushi_apps = runnable_application(@data_set.headers, refresh)
      @sushi_apps_category = sushi_apps.map{|app| app.analysis_category}.uniq.sort
      @sushi_apps = {}
      sushi_apps.sort_by{|app| app.class_name.to_s}.each do |app|
        @sushi_apps[app.analysis_category] ||= []
        @sushi_apps[app.analysis_category] << app.class_name.to_s
      end
      @data_set.runnable_apps = @sushi_apps
      @data_set.save
    end
  end
  def show
    view_context.project_init
    @fgcz = SushiFabric::Application.config.fgcz?
    # switch project (from job_monitoring)
    if project = params[:project]
      session[:project] = project.to_i
    end
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

    @data_set = DataSet.find_by_id(params[:id])

    # check some properties
    if session[:employee]
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
              file_path = File.join(SushiFabric::GSTORE_DIR, file)
              @sample_path << File.dirname(file)
              @file_exist[file] = File.exist?(file_path)
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

    # update num_samples
    if @data_set.num_samples.to_i != sample_count
      @data_set.num_samples = sample_count
    end

    if !@data_set.refreshed_apps and @data_set.runnable_apps.empty?
      @data_set.refreshed_apps = true
      @data_set.save
      set_runnable_apps(false)
    end
    if @file_exist.values.inject{|a,b| a and b}
      @sushi_apps = @data_set.runnable_apps
      @sushi_apps_category = @sushi_apps.keys.sort
    end
  end
  def refresh_apps
    set_runnable_apps
    show
  end
  def edit
    show
  end
  def trace_treeviews(root, data_set, parent_id, project_number)
    data_set_id = data_set.id
    node_text = if root.empty?
             "<b>" + data_set.data_sets.length.to_s+" "+data_set.name+"</b> "+" <small><font color='gray'>"+data_set.comment.to_s+"</font></small>"
           else
              data_set.data_sets.length.to_s+" "+data_set.name+" "+" <small><font color='gray'>"+data_set.comment.to_s+"</font></small>"
           end
    node = {"id" => data_set_id, 
            "text" => node_text,
            "parent" => parent_id,
            "state" => {"opened":true},
            "a_attr" => {"href"=>"/data_set/p#{project_number}/#{data_set_id}", 
                         "onclick"=>"window.open('/data_set/p#{project_number}/#{data_set_id}')"}
            }
    root << node
    data_set.data_sets.each do |child|
      if child.project.number==project_number
        trace_treeviews(root, child, data_set.id, project_number)
      end
    end
  end
  def partial_treeviews
    root = []
    if top_data_set_id = params[:format]
      data_set = DataSet.find_by_id(top_data_set_id.to_i)
      if children = data_set.data_sets and children.length > 0
        trace_treeviews(root, data_set, "#", data_set.project.number)
      end
    end
    render :json => root.reverse
  end
  def whole_treeviews
    @project = Project.find_by_number(session[:project].to_i)
    root = []
    project_dataset_ids = Hash[*(@project.data_sets.map{|data_set| [data_set.id, true]}.flatten)]
    @project.data_sets.each do |data_set|
      node = {"id" => data_set.id, 
              "text" => data_set.data_sets.length.to_s+" "+data_set.name+" <small><font color='gray'>"+data_set.comment.to_s+"</font></small>",
              "a_attr" => {"href"=>"/data_set/p#{@project.number}/#{data_set.id}", 
                           "onclick"=>"window.open('/data_set/p#{@project.number}/#{data_set.id}')"}
              }
      if parent = data_set.data_set and project_dataset_ids[parent.id]
        node["parent"] = parent.id
      else
        node["parent"] = "#"
      end
      root << node
    end
    
    render :json => root.reverse
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
              @data_set_id = save_data_set(data_set, headers, rows, current_user)
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
        data_set_tsv.each do |row|
          unless row.fields.join.strip.empty?
            rows << row.fields
          end
        end
        unless headers.include?(nil)
          @data_set_id = save_data_set(data_set, headers, rows, current_user)
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

        unless session[:off_bfabric_registration]
          data_set = DataSet.find_by_id(@data_set_id)
          pid = Process.fork do
            Process.fork do
              data_set.register_bfabric
            end # grand-child process
          end # child process
          Process.waitpid pid
        end
      end
    end

    redirect_to :controller => "data_set"
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
        multi_data_sets = false
        open(tsv.path) do |input|
          while line=input.gets
            if line =~ /ProjectNumber/
              multi_data_sets = true
              break
            end
          end
        end
        
        if multi_data_sets
          csv = CSV.readlines(tsv.path, :col_sep=>"\t")
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
                @data_set_id = save_data_set(data_set, headers, rows, current_user)
              else
                @warning = 'There must be a blank column. Please check it. Import is incomplete.'
              end
            end
            if row.empty?
              data_set = []
              headers = []
              rows = []
            end
          end
        else
          data_set_tsv = CSV.readlines(tsv.path, :headers => true, :col_sep=>"\t")

          data_set = []
          headers = data_set_tsv.headers
          rows = []
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
            end
          end
            unless headers.include?(nil)
              @data_set_id = save_data_set(data_set, headers, rows, current_user)
            else
              @warning = 'There must be a blank column. Please check it. Import is incomplete.'
            end
        end
      end

      if @data_set_id
        set_runnable_apps(false)

        unless session[:off_bfabric_registration]
          data_set = DataSet.find_by_id(@data_set_id)
          pid = Process.fork do
            Process.fork do
              data_set.register_bfabric
            end # grand-child process
          end # child process
          Process.waitpid pid
        end
      elsif file = params[:file] and tsv = file[:name]
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
      @option = params[:option]

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
        @command = @@workflow_manager.delete_command(target)
        if @option[:delete] == 'only_gstore' or @option[:delete] == 'also_gstore'
          @command_log = `#{@command}`
          if request = @command_log.split and request_no = request[4]
            @greq_status_command = "g-req status #{request_no}"
          end
        end
      end

      # delete data in sushi
      if @option[:delete] == 'only_sushi' or @option[:delete] == 'also_gstore'
        @data_set.samples.each do |sample|
          sample.delete
        end
        @deleted_data_set = @data_set.delete
      else
        @deleted_data_set = @data_set
      end
    end
  end
  def multi_destroy
    @data_set_ids = params[:option][:data_set_ids].split(',')
    @commands = []
    @command_logs = []
    @deleted_data_sets = []
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
  end
  def job_parameter
    @data_set = if id = params[:format]
                  DataSet.find_by_id(id)
                end
    @sample_path = if @data_set 
                     sample_path(@data_set)
                   end
    @parameters_tsv = if @sample_path
                        File.join(SushiFabric::GSTORE_DIR, @sample_path, 'parameters.tsv')
                      end
    @parameters = {}
    if @parameters_tsv and File.exist?(@parameters_tsv)
      File.readlines(@parameters_tsv).each do |line|
        header, *values = line.chomp.split
        @parameters[header] = values.join(" ")
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
    @command = @@workflow_manager.delete_command(target)
    @command_log = `#{@command}`
  end
  def register_bfabric
    data_set = DataSet.find_by_id(params[:id])
    pid = Process.fork do
      Process.fork do
        data_set.register_bfabric("only_one")
      end # grand-child process
    end # child process
    Process.waitpid pid
    redirect_to :controller => "data_set", :action => "show"
  end
end
