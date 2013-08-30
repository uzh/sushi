class SushiApplicationController < ApplicationController
  def index
    @sushi_apps = all_sushi_applications.map{|app| eval(app.gsub(/\.rb/,'').gsub(/\.sh/,'')).new}
    @sushi_apps.each do |app|
      app.instance_variable_set(:@dataset, {})
    end
    @sushi_apps = @sushi_apps.sort_by{|app| app.analysis_category}
  end
end
