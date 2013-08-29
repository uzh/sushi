class SampleController < ApplicationController
  def show
    @data_set = DataSet.find_by_id(params[:id])

    # check real data
    @file_exist = {}
    @sample_path = []
    @data_set.samples.each do |sample|
      sample.to_hash.each do |header, file|
        if header =~ /\[File\]/
          file_path = File.join(GSTORE_DIR, file)
          @sample_path << File.dirname(file)
          @file_exist[file] = File.exist?(file_path)
        else
          @file_exist[file] = true
        end
      end
    end
    @sample_path.uniq!

    if @file_exist.values.inject{|a,b| a and b}
      @sushi_apps = runnable_application(@data_set.headers)
    end
  end
end
