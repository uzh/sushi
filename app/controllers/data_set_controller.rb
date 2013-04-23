class DataSetController < ApplicationController
  def index
    params[:project] = session[:project]

    if session[:project] 
      unless @project = Project.find_by_number(session[:project].to_i)
        @project = Project.new
        @project.number = session[:project].to_i
        @project.save
      end
      if file = params[:file] and csv = file[:name]
        data_set_csv = CSV.readlines(csv.path, :headers => true)

        @data_set = DataSet.new
        if dataset = params[:dataset] and dataset_name = dataset[:name]
          @data_set.name = dataset_name
        else
          @data_set.name = 'DataSet ' + (DataSet.all.length+1).to_s
        end
        data_set_csv.each do |row|
          sample = Sample.new
          sample.key_value = row.to_hash.to_s
          sample.save unless sample.saved?
          @data_set.samples << sample
        end
        @data_set.md5 = @data_set.md5hexdigest
        unless @data_set.saved?
          @project.data_sets << @data_set
          @data_set.save 
        end
      end
    end
  end

  def save_as_csv
    project_dir = Dir.pwd
    data_set_csv = File.join(project_dir, 'public/test_dataset.csv')
    if id = params[:id] and data_set = DataSet.find_by_id(id)
=begin
      CSV.open(data_set_csv, 'w') do |csv|
        csv << data_set.headers
        new_data_set.new_samples.each do |sample|
          row = []
          row_hash = sample.to_hash
          data_set.headers.each do |header|
            row << row_hash[header]
          end 
          csv << row
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
