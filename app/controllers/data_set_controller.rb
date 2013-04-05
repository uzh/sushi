class DataSetController < ApplicationController
  def index
    if file = params[:csv_file] and file_io = file[:name]
      csv = CSV.readlines(file_io.path, :headers => true)
      @data_set = DataSet.new
      csv.each do |row|
        sample = row.to_hash.to_s
        new_sample = Sample.new
        new_sample.key_value = row.to_hash.to_s
        #new_sample.save unless new_sample.saved?
        @data_set.samples << new_sample
      end
      if data_set = params[:data_set] and data_set_name = data_set[:name]
        @data_set.name = data_set_name
      end
      @data_set.save unless @data_set.saved?
    end
    @data_sets = DataSet.all
    @samples = Sample.all
  end
  def treeviews
    tree = []
    node_list = {}
    DataSet.all.each do |data_set|
      node = {"id" => data_set.id, "text" => data_set.id.to_s+' '+data_set.name, 'path' => '', "expanded" => true, "classes" => 'file', "hasChildren" => false, "children" => []}
      node_list[data_set.id] = node
      if parent = data_set.data_set
        node_list[parent.id]['children'] << node
      else
        tree << node
      end
    end
    render :json => tree
  end
  def edit
=begin
    if delete_parent_id = params[:delete_parent_id]
      parent_id = delete_parent_id.to_i
      parent = DataSet.find_by_id(parent_id)
      child_id =  params[:id].to_i
      child = DataSet.find_by_id(child_id)
      parent.data_sets.delete(child)
    end
=end
    @data_sets = DataSet.all
    @data_set = DataSet.find_by_id(params[:id])
    @samples = Sample.all
    @params = params
  end
  def add
    @data_sets = DataSet.all
    @samples = Sample.all
    @params = params
    render 'data_set/edit'
  end
  def delete
    if sample_ids = params[:sample_ids]
      sample_ids.each do |id|
        Sample.find(id.to_i).destroy
      end
    end
    if data_set_ids = params[:data_set_ids]
      data_set_ids.each do |id|
        DataSet.find(id.to_i).destroy
      end
    end
    @data_sets = DataSet.all
    @samples = Sample.all
    @params = params
  end
  def add_or_delete
    params[:method] = params[:method].strip
    data_set_id = nil
    sample_id = nil
    if params[:method] == 'add'
      params[:samples].each do |sample|
        id = sample.split(/,/)
        data_set_id = id[0].to_i
        sample_id = id[1].to_i
      end
    else
      params[:data_sets].each do |data_set|
        id = data_set.split(/,/)
        data_set_id = id[0].to_i
        sample_id = id[1].to_i
      end
    end
    @data_sets = DataSet.all
    @samples = Sample.all
    @params = params
    render 'data_set/edit'
  end
  
  def create
    if data_set_name = params[:data_set_name] and sample_resource_ids = params[:sample_resource_ids].map{|i| i.to_s}
      samples = session[:basket].select{|sample| sample_resource_ids.include?(sample.resource_id)}
      samples.each_with_index do |sample,i|
        unless sample.id
          sample.name = params[:sample_names][i]
          sample.save
        end
      end
      ds = DataSet.new
      ds.name = data_set_name
      ds.save
      if data_set = params[:data_set] and parent = data_set[:parent] and !parent.empty?
        id = parent.split(/ /)[0].to_i
        parent_data_set = DataSet.find_by_id(id)
        parent_data_set.data_sets << ds 
      end
    end

    @sushis = session[:basket]

    @data_sets = DataSet.all
  end
end
