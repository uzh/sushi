class DataListController < ApplicationController
  def index
    @data_lists = DataList.all
  end
  def add
    data_list = DataList.new
    data_list.data_set_id = 2
    data_list.sample_id = 1
    data_list.save
    render :text => 'A DataList is added'
  end
  def delete
    id = DataList.all.last.id
    DataList.delete(id)
    render :text => "DataList (id=#{id}) is deleted"
  end
end
