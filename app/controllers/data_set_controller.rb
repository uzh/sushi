class DataSetController < ApplicationController
  def index
    @data_sets = DataSet.all
    @data_lists = DataList.all
    @samples = Sample.all
  end
  def edit
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
    id = DataList.all.last.id
    DataList.delete(id)
    @data_sets = DataSet.all
    @data_lists = DataList.all
    @samples = Sample.all
    @params = params
    render 'data_set/edit'
  end
end
