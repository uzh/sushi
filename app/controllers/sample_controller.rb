class SampleController < ApplicationController
  def show
    @data_set = DataSet.find_by_id(params[:id])

    # update column names
    # assuming values are not different, 
    # in other words, this should be done after editing values
    new_headers = params[:sample_headers]
    current_headers = Hash[*@data_set.samples.first.to_hash.keys.map{|header| [header.split.first, header]}.flatten]
    if new_headers and new_headers!=current_headers
      @data_set.samples.each_with_index do |sample, i|
        new_sample = {}
        sample.to_hash.each do |header, value|
          header_without_tag = header.split.first
          new_header = new_headers[header_without_tag]
          new_sample[new_header] = value
        end
        @data_set.samples[i].key_value = new_sample.to_s
        @data_set.samples[i].save
      end
      @data_set.md5 = @data_set.md5hexdigest
      @data_set.save
    end

    # update values
=begin
    params.select{|key| key=~/sample_\d+/}.sort_by{|key| key.match(/sample_(\d+)/)[1].to_i}.each do |sample_no|
      if sample = params[sample_no] and i = sample_no.split(/_/)[1].to_i
        sample = eval(sample)
        new_sample = {}
        @data_set.samples[i].to_hash.keys.each do |header|
          if new_value = sample[saheader.split.first]
            new_sample[header] = new_value
          end
        end
      end
    end
=end
  end
  def edit
    @data_set = DataSet.find_by_id(params[:id])
  end
end
