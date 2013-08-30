class SampleController < ApplicationController
  def show
    @data_set = DataSet.find_by_id(params[:id])

    # update values
    @data_set.samples.each_with_index do |sample, i|
      current_sample = Hash[*sample.to_hash.map{|key, value| [key.split.first,value]}.flatten]
      if edit_sample = params["sample_#{i}"] and current_sample.to_s!=edit_sample.to_s
        new_sample = {}
        sample.to_hash.each do |header, value|
          header_without_tag = header.to_s.split.first 
          new_sample[header] = edit_sample[header_without_tag]
        end
        @data_set.samples[i].key_value = new_sample.to_s
        @data_set.samples[i].save
        @data_set.md5 = @data_set.md5hexdigest
        @data_set.save
      end
    end

    # add new row
    current_headers = Hash[*@data_set.samples.first.to_hash.keys.map{|header| [header.to_s.split.first, header]}.flatten]
    if add_sample = params[:sample_new]
      new_sample = {}
      add_sample.each do |key, value|
        header = current_headers[key]
        if header =~ /\[File\]/ and sample = @data_set.samples.first and sample = sample.to_hash[header]
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
    end

    # del row
    if edit_option = params[:edit_option] and row = edit_option[:del_row]
      @data_set.samples[row.to_i].delete
      @data_set.md5 = @data_set.md5hexdigest
      @data_set.save
      @data_set = DataSet.find_by_id(params[:id])
    end

    # update column names
    # assuming values are not different, 
    # in other words, this should be done after editing values
    new_headers = params[:sample_headers]
    #current_headers = Hash[*@data_set.samples.first.to_hash.keys.map{|header| [header.to_s.split.first, header]}.flatten]
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
    end
  end
  def edit
    @data_set = DataSet.find_by_id(params[:id])
    @edit_option = params[:edit_option]
  end
end
