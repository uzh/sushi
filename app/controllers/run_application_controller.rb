class RunApplicationController < ApplicationController
  def index
    @job_scripts = Dir['lib/*.rb'].sort.select{|script| script !~ /sushiApp\.rb/}.to_a.map{|script| script.gsub(/\.rb/,'')}
    if project_number = session[:project] and project = Project.find_by_number(project_number.to_i)
      @data_sets = project.data_sets
    end
  end
  def set_parameters
  end
end
