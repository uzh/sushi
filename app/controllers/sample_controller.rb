class SampleController < ApplicationController
  def show
    @data_set = DataSet.find_by_id(params[:id])

    # update values
    @hoge=false
    @data_set.samples.each_with_index do |sample, i|
      current_sample = Hash[*sample.to_hash.map{|key, value| [key.split.first,value]}.flatten]
      if edit_sample = params["sample_#{i}"] and current_sample.to_s!=edit_sample.to_s
        @hoge=true
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

    # update column names
    # assuming values are not different, 
    # in other words, this should be done after editing values
    new_headers = params[:sample_headers]
    current_headers = Hash[*@data_set.samples.first.to_hash.keys.map{|header| [header.to_s.split.first, header]}.flatten]
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
  end
  def edit
    @data_set = DataSet.find_by_id(params[:id])
  end
end
