class SampleController < ApplicationController
  include SushiFabric
  def show
    @data_set = DataSet.find_by_id(params[:id])

    # update values
    samples = []
    @data_set.samples.each_with_index do |sample, i|
      current_sample = Hash[*sample.to_hash.map{|key, value| [key.split.first,value]}.flatten]
      if edit_sample = params["sample_#{i}"] and current_sample.to_s!=edit_sample.to_s
        new_sample = {}
        sample.to_hash.each do |header, value|
          header_without_tag = header.to_s.split.first 
          new_sample[header] = edit_sample[header_without_tag]
        end
        if params[:edit]
          @data_set.samples[i].key_value = new_sample.to_s
          @data_set.samples[i].save
          @data_set.md5 = @data_set.md5hexdigest
          @data_set.save
        elsif params[:edit_save_as_child]
          samples << new_sample.values
        end
      elsif params[:edit_save_as_child]
        samples << sample.to_hash.values
      end
    end
    if params[:edit]
      save_dataset_tsv_in_gstore(@data_set)
    end

    if params[:edit_save_as_child] and del_rows = session[:del_rows]
      del_rows.sort.reverse.each do |i|
        samples.delete_at(i)
      end
    end

    # save as a child dataset
    if params[:edit_save_as_child]    
      data_set = []
      data_set << "DataSetName"
      data_set << if child_dataset = params[:child_dataset] and child_name = child_dataset[:name]
                    child_name
                  else
                    "Child_#{@data_set.name}"
                  end
      project = Project.find_by_number(session[:project].to_i)
      data_set << "ProjectNumber" << project.number
      data_set << "ParentID" << @data_set.id 
      headers = @data_set.headers
      if data_set_id = save_data_set(data_set, headers, samples, current_user)
        @data_set = DataSet.find_by_id(data_set_id)
        @data_set.child = true
        @data_set.save
        save_dataset_tsv_in_gstore(@data_set)
      end
    end

    # add new row
    current_headers = Hash[*@data_set.samples.first.to_hash.keys.map{|header| [header.to_s.split.first, header]}.flatten]
    if add_sample = params[:sample_new]
      new_sample = {}
      add_sample.each do |key, value|
        header = current_headers[key]
        if header.tag?('File') and sample = @data_set.samples.first and sample = sample.to_hash[header]
          new_sample[header] = File.join(File.dirname(sample), value)
        else
          new_sample[header] = value
        end
      end
      sample = Sample.new
      sample.key_value = new_sample.to_s
      sample.save
      @data_set.samples << sample
      @data_set.md5 = @data_set.md5hexdigest
      @data_set.save
      save_dataset_tsv_in_gstore(@data_set)
    end

    # del row
    if params[:edit_save] and del_rows = session[:del_rows]
      del_rows.sort.reverse.each do |i|
        @data_set.samples.delete(@data_set.samples[i])
      end
      @data_set.num_samples = @data_set.samples.length - del_rows.length
      @data_set.md5 = @data_set.md5hexdigest
      @data_set.save
      save_dataset_tsv_in_gstore(@data_set)
    end

    # update column names
    # assuming values are not different, 
    # in other words, this should be done after editing values
    new_headers = params[:sample_headers]
    if new_headers and new_headers!=current_headers
      @data_set.samples.each_with_index do |sample, i|
        new_sample = {}
        sample.to_hash.each do |header, value|
          header_without_tag = header.split.first 
          if new_header = new_headers[header_without_tag] and !new_header.to_s.empty?
            new_sample[new_header] = value
          else
            new_sample[header] = value
          end
        end
        @data_set.samples[i].key_value = new_sample.to_s
        @data_set.samples[i].save
      end
      @data_set.md5 = @data_set.md5hexdigest
      @data_set.save
      save_dataset_tsv_in_gstore(@data_set)
    end

    # add new column
    if new_header = params[:new_header] and new_header_name = new_header[:name]
      @data_set.samples.each_with_index do |sample, i|
        new_sample = sample.to_hash
        new_value = params[:new_col][i.to_s]
        new_sample[new_header_name]=new_value
        @data_set.samples[i].key_value = new_sample.to_s
        @data_set.samples[i].save
      end
      @data_set.md5 = @data_set.md5hexdigest
      @data_set.save
      save_dataset_tsv_in_gstore(@data_set)
    end

    # delete column
    if edit_option = params[:edit_option] and header = edit_option[:del_col]
      @data_set.samples.each_with_index do |sample, i|
        new_sample = sample.to_hash
        new_sample.delete(header)
        @data_set.samples[i].key_value = new_sample.to_s
        @data_set.samples[i].save
      end
      @data_set.md5 = @data_set.md5hexdigest
      @data_set.save
      save_dataset_tsv_in_gstore(@data_set)
    end

    # delete rows session clear
    if session[:del_rows]
      session[:del_rows] = nil
    end
  end
  def edit
    @data_set = DataSet.find_by_id(params[:id])
    @edit_option = params[:edit_option]
  end
  def multiedit
    @data_set = DataSet.find_by_id(params[:id])
    @edit_option = params[:edit_option]
    if del_row = @edit_option[:del_row]
      session[:del_rows] ||= []
      session[:del_rows] << del_row.to_i
      session[:del_rows].uniq!
    end
  end
end
