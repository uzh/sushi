class SushiApplicationController < ApplicationController
  def index
    #@sushi_apps = all_sushi_applications.map{|app| eval(app.gsub(/\.rb/,'').gsub(/\.sh/,'')).new}
    #@sushi_apps = all_sushi_applications.map{|app| ''}
    #@sushi_apps = []
    #@sushi_apps.each do |app|
    #  app.instance_variable_set(:@dataset, {})
    #end
    #@sushi_apps = @sushi_apps.sort_by{|app| app.analysis_category}

    @sushi_apps = SushiApplication.all.sort_by{|app| app.analysis_category}
  end
  def refresh
    sushi_apps = all_sushi_applications
    sushi_apps.each do |app|
      class_name = app.gsub(/\.rb/,'').gsub(/\.sh/,'')
      unless SushiApplication.find_by_class_name(class_name)
        sushi_app_instance = eval(class_name).new
        sushi_app_instance.instance_variable_set(:@dataset, {})
        sushi_app_instance.instance_variable_set(:@result_dir, '')
        sushi_app_entry = SushiApplication.new
        sushi_app_entry.class_name = class_name
        sushi_app_entry.analysis_category = sushi_app_instance.analysis_category
        sushi_app_entry.required_columns = sushi_app_instance.required_columns
        sushi_app_entry.next_dataset_keys = sushi_app_instance.next_dataset.keys
        sushi_app_entry.save
      end
    end
    @sushi_apps = SushiApplication.all.sort_by{|app| app.analysis_category}
    render action: "index"
  end
end
