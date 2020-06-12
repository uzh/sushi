require 'json'
require 'csv'

class ApiController < ApplicationController
  def save_dataset project_id, dataset_name, headers, rows
    project = Project.find_by_number project_id
    if not project
      project = Project.new
      project.number = project_id
      project.save
    end 
    
    data_set = DataSet.new
    data_set.name = dataset_name
    data_set.project = project
    sample_hash = {}
    rows.each do |row|
      headers.each_with_index do |header, i|
        sample_hash[header] = row[i]
      end
      sample = Sample.new
      sample.key_value = sample_hash.to_s
      sample.save unless sample.saved?
      data_set.samples << sample
    end
    
    data_set.md5 = data_set.md5hexdigest
    data_set.save
     
    return data_set.id.to_i
  end
  
  def index
    if not request.method.eql? "POST"
      render :plain => "Please POST data in JSON format.", :status => :bad_request
      return
    end
     
    begin
      json = JSON.parse request.body.read
    rescue
      render :plain => "Error reading JSON data.", :status => :bad_request
      return
    end
    
    if not json.has_key? "project"
      render :plain => "Please provide a project number: 'project'", :status => :bad_request
      return
    end
    
    project_id = json["project"].to_i
    if project_id == 0
      render :plain => "Bad project number!", :status => :bad_request
      return
    end
    
    if not json.has_key? "name"
      render :plain => "Please provide a dataset name: 'name'", :status => :bad_request
      return
    end
    dataset_name = json["name"]
    
    if not json.has_key? "path"
      render :plain => "Please provide a dataset file: 'path'", :status => :bad_request
      return
    end
    
    dataset_file_path = json["path"]
    rows = []
    begin
      rows = CSV.readlines dataset_file_path, :col_sep=>"\t"
    rescue
      render :plain => "Cannot read the dataset file: "+dataset_file_path, :status => :bad_request
      return
    end
    
    if rows.length < 2
      render :plain => "Dataset must contain one header line and at least one sample line.", :status => :bad_request
      return
    end
    
    headers = rows.shift
    ds_id = save_dataset project_id, dataset_name, headers, rows
    if ds_id == 0
      render :plain => "Cannot save the dataset!", :status => :internal_server_error
      return
    end
    
    response = {}
    response['message'] = "OK"
    response['data_set_id'] = ds_id
    render :json => response, :status => :ok
  end
end
