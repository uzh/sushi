class RunApplicationController < ApplicationController
  def index
    @job_scripts = Dir['public/*.sh'].sort.to_a
    @data_sets = DataSet.all
  end
end
