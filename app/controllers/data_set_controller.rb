class DataSetController < ApplicationController
  def index
    @data_sets = DataSet.all
    @data_lists = DataList.all
    @samples = Sample.all
  end
  def edit
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
    data_list = DataList.find_by_data_set_id_and_sample_id(params[:data_set_id].to_i,params[:sample_id].to_i)
    DataList.delete(data_list.id)
    @data_sets = DataSet.all
    @data_lists = DataList.all
    @samples = Sample.all
    @params = params
    render 'data_set/edit'
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
    if params.include? :data_set_name
      if params.include? :sample_ids
        ds = DataSet.new
        ds.note = params[:data_set_name]
        ds.save
        
        params[:sample_ids].each do |s_id|
          dl = DataList.new
          dl.data_set_id = ds.id
          dl.sample_id = s_id.to_i
          dl.save
        end
      end
    end

    @sushis = Array.new
    if session.include? :basket
      session[:basket].each do |r_id|
        @sushis << Sample.find_by_resource_id(r_id)
      end 
    end

    @data_sets = DataSet.all
  end
end
