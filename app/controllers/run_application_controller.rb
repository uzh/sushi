class RunApplicationController < ApplicationController
  def index
    @job_scripts = Dir['lib/*.rb'].sort.select{|script| script !~ /sushiApp\.rb/}.to_a.map{|script| script.gsub(/\.rb/,'')}
    @data_sets = DataSet.all
  end
  def set_parameter
  end
end
