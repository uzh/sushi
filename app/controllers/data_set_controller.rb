class DataSetController < ApplicationController

=begin
 def save_data_set(data_set_arr, headers, rows)
    data_set_hash = Hash[*data_set_arr]
    if project = Project.find_by_number(data_set_hash['ProjectNumber'].to_i)
      data_set = DataSet.new
      data_set.name = data_set_hash['Name']
      data_set.project = project
      if parent_id = data_set_hash['ParentID'] and parent_data_set = DataSet.find_by_id(parent_id.to_i)
        data_set.data_set = parent_data_set
      end

      sample_hash = {}
      rows.each do |row|
        headers.each_with_index do |header, i|
         sample_hash[header]=row[i]
        end
        sample = Sample.new
        sample.key_value = sample_hash.to_s
        unless sample.saved?
          sample.save unless sample.saved?
        end
        data_set.samples << sample
      end

      data_set.md5 = data_set.md5hexdigest
      unless data_set.saved?
        project.data_sets << data_set
        parent_data_set.data_sets << data_set if parent_data_set
        data_set.save
      end

    end
  end
=end
  def index
    params[:project] = session[:project]

    if session[:project] 
      unless @project = Project.find_by_number(session[:project].to_i)
        @project = Project.new
        @project.number = session[:project].to_i
        @project.save
      end

      if file = params[:file] and tsv = file[:name]
        data_set_set = false
        open(tsv.path) do |input|
          while line=input.gets
            if line =~ /ProjectNumber/
              data_set_set = true
              break
            end
          end
        end
        
        if data_set_set

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

          @data_set = DataSet.new
          @data_set.project = @project
          if dataset = params[:dataset] and dataset_name = dataset[:name]
            @data_set.name = dataset_name
          else
            @data_set.name = 'DataSet ' + (DataSet.all.length+1).to_s
          end
          data_set_tsv.each do |row|
            sample = Sample.new
            sample.key_value = row.to_hash.to_s
            sample.save unless sample.saved?
            @data_set.samples << sample
          end

          if parent = params[:parent] and parent_id = parent[:id] and @parent_data_set = DataSet.find_by_id(parent_id.to_i)
            @data_set.data_set = @parent_data_set
          end

          @data_set.md5 = @data_set.md5hexdigest
    
          unless @data_set.saved?
            @project.data_sets << @data_set
            @parent_data_set.data_sets << @data_set if @parent_data_set
            @data_set.save 
          end
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
