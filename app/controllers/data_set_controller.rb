class DataSetController < ApplicationController
  def index
    @data_sets = DataSet.all
    @data_lists = DataList.all
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
    if delete_parent_id = params[:delete_parent_id]
      parent_id = delete_parent_id.to_i
      parent = DataSet.find_by_id(parent_id)
      child_id =  params[:id].to_i
      child = DataSet.find_by_id(child_id)
      parent.data_sets.delete(child)
    end

    @data_sets = DataSet.all
    @data_set = DataSet.find_by_id(params[:id])
    @data_lists = DataList.all
    @samples = Sample.all
    @params = params
  end
  def add
    data_list = DataList.new
    data_list.data_set_id = params[:data_set_id].to_i
    data_list.sample_id = params[:sample_id].to_i
    data_list.save
    @data_sets = DataSet.all
    @data_lists = DataList.all
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
    if data_list_ids = params[:data_list_ids]
      data_list_ids.each do |id|
        DataList.find(id.to_i).destroy
      end
    end
    if data_set_ids = params[:data_set_ids]
      data_set_ids.each do |id|
        DataSet.find(id.to_i).destroy
      end
    end
    @data_sets = DataSet.all
    @data_lists = DataList.all
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
        data_list = DataList.new
        data_list.data_set_id = data_set_id
        data_list.sample_id = sample_id
        data_list.save
      end
    else
      params[:data_sets].each do |data_set|
        id = data_set.split(/,/)
        data_set_id = id[0].to_i
        sample_id = id[1].to_i
        data_list = DataList.find_by_data_set_id_and_sample_id(data_set_id,sample_id)
      DataList.delete(data_list.id)
      end
    end
    @data_sets = DataSet.all
    @data_lists = DataList.all
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
     
      samples.each do |sample|
        dl = DataList.new
        dl.data_set_id = ds.id
        dl.sample_id = sample.id
        dl.save
      end
    end

    @sushis = session[:basket]

    @data_sets = DataSet.all
  end
end
