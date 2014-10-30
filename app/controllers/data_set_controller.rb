class DataSetController < ApplicationController
  include SushiFabric
  def index
    @project = Project.find_by_number(session[:project].to_i)

    @sample_available = {}
    if @project
      @project.data_sets.each do |data_set|
        data_set.samples.each do |sample|
          flag = true
          sample.to_hash.each do |header, file|
            if header and header.tag?('File') and file and file_path = File.join(SushiFabric::GSTORE_DIR, file) and !File.exist?(file_path)
              flag = false
              break
            end
          end
          @sample_available[data_set] ||= 0
          @sample_available[data_set] += 1 if flag
        end
      end
    end
  end
  def script_log
    @data_set = if id = params[:format]
                  DataSet.find_by_id(id)
                end
  end
  def show
    @fgcz = (`hostname`.chomp =~ /fgcz-s-034/)
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

    # check real data
    @file_exist = {}
    @sample_path = []
    @sample_invalid_name = {}
    if @data_set
      @data_set.samples.each do |sample|
        sample.to_hash.each do |header, file|
          if header and (header.tag?('File') or header.tag?('Link')) 
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

    if @file_exist.values.inject{|a,b| a and b}
      sushi_apps = runnable_application(@data_set.headers)
      sushi_apps = sushi_apps.map{|app| eval(app.gsub(/\.rb/,'').gsub(/\.sh/,'')).new}
      @sushi_apps_category = sushi_apps.map{|app| app.analysis_category}.uniq.sort
      @sushi_apps = {}
      sushi_apps.sort_by{|app| app.class.to_s}.each do |app|
        @sushi_apps[app.analysis_category] ||= []
        @sushi_apps[app.analysis_category] << app.class.to_s
      end
    end
  end
  def edit
    show
  end
  def treeviews
    @project = Project.find_by_number(session[:project].to_i)
    tree = []
    node_list = {}
    root = []
    top_nodes = []
    @project.data_sets.each do |data_set|
      node = {"id" => data_set.id, 
              "text" => data_set.data_sets.length.to_s+
              " <a href='/data_set/#{data_set.id}'>"+data_set.name+'</a>'+
              ' <span style="color:gray;font-size:smaller">['+data_set.headers.join(',')+']</span>', 
              'path' => '', 
              "expanded" => false, 
              "classes" => 'file', 
              "hasChildren" => false, 
              "children" => []}
      node_list[data_set.id] = node
      if parent = data_set.data_set
        node_list[parent.id]['children'] << node
      else
        top_nodes << node
      end
      if data_set.id == params[:format].to_i
        root << node
      end
    end
    root = top_nodes.reverse if root.empty?
    tree.concat root
    render :json => tree
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
          elsif !row.empty?
            rows << row
          else
            save_data_set(data_set, headers, rows)
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
          rows << row.fields
        end
        save_data_set(data_set, headers, rows)
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
            elsif !row.empty?
              rows << row
            else
              save_data_set(data_set, headers, rows)
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
            rows << row.fields
          end
          save_data_set(data_set, headers, rows)
        end
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

  def destroy
    @fgcz = (`hostname` =~ /fgcz-s-034/)
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
        @command = `wfm_delete_command -f #{File.join(SushiFabric::GSTORE_DIR, @sample_path.first)} -d #{SushiFabric::WORKFLOW_MANAGER}`
        if @option[:delete] == 'also_gstore'
          @command_log = `#{@command}`
          if request = @command_log.split and request_no = request[4]
            @greq_status_command = "g-req status #{request_no}"
          end
        end
      end

      # delete data in sushi
      @data_set.samples.each do |sample|
        sample.delete
      end
      @deleted_data_set = @data_set.delete

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
end
