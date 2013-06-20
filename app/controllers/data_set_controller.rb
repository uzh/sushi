class DataSetController < ApplicationController
  def index
    @project = Project.find_by_number(session[:project].to_i)
    if new_data_set = params[:data_set] and name = new_data_set[:name] and id = new_data_set[:id]
      data_set = DataSet.find_by_id(id)
      data_set.name = name
      data_set.save
    end
  end
  def show
    @data_set = DataSet.find_by_id(params[:id])
  end
  def edit
    @project = Project.find_by_number(session[:project].to_i)
  end
  def treeviews
    @project = Project.find_by_number(session[:project].to_i)
    tree = []
    node_list = {}
    @project.data_sets.each do |data_set|
      node = {"id" => data_set.id, "text" => data_set.id.to_s+" <a href='/data_set/#{data_set.id}'>"+data_set.name+'</a>', 'path' => '', "expanded" => true, "classes" => 'file', "hasChildren" => false, "children" => []}
      node_list[data_set.id] = node
      if parent = data_set.data_set
        node_list[parent.id]['children'] << node
      else
        tree << node
      end
    end
    render :json => tree
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
    project_dir = Dir.pwd
    data_set_tsv = File.join(project_dir, 'public/test_dataset.tsv')
    if id = params[:id] and data_set = DataSet.find_by_id(id)
=begin
      CSV.open(data_set_tsv, 'w', :col_sep=>"\t") do |tsv|
        tsv << data_set.headers
        new_data_set.new_samples.each do |sample|
          row = []
          row_hash = sample.to_hash
          data_set.headers.each do |header|
            row << row_hash[header]
          end 
          tsv << row
        end
      end
=end
    end

    file = params[:file][:name]
    file_name = file.original_filename
    file_size = file.size
    #path = File.join('public',file_name) 
    #UploadFile.create(:filename => file_name, :filesize => file_size, :filepath => path)
    #File.open(path, "w") { |f| f.write(file.read) }
    render :text => "params: #{params} file.class: #{file.path} file_path: #{file_name}"
  end
end
