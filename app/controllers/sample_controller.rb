class SampleController < ApplicationController
  def show
    @data_set = DataSet.find_by_id(params[:id])

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
