class SushiApplicationController < ApplicationController
  def index
    @sushi_apps = SushiApplication.all.sort_by{|app| app.analysis_category}
  end
  def refresh
    refresh_sushi_application
    @sushi_apps = SushiApplication.all.sort_by{|app| app.analysis_category}
    render action: "index"
  end
end
