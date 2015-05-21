class SushiApplicationController < ApplicationController
  def index
    @sushi_apps = SushiApplication.all.sort_by{|app| app.analysis_category}
  end
  def refresh

    # new check
    sushi_apps = all_sushi_applications
    lib_dir = File.expand_path('../../../lib', __FILE__)
    sushi_apps.select{|app| app =~ /\.rb$/}.each do |app|
      class_name = app.gsub(/\.rb/,'')
      updated_at = File.stat(File.join(lib_dir, app)).mtime
      unless sushi_app = SushiApplication.find_by_class_name(class_name) and updated_at < sushi_app.updated_at
        if sushi_app 
          load File.join(lib_dir, app)
        end
        sushi_app_instance = eval(class_name).new
        sushi_app_instance.instance_variable_set(:@dataset, {})
        sushi_app_instance.instance_variable_set(:@result_dir, '')
        sushi_app_entry = (sushi_app || SushiApplication.new)
        sushi_app_entry.class_name = class_name
        sushi_app_entry.analysis_category = sushi_app_instance.analysis_category
        sushi_app_entry.required_columns = sushi_app_instance.required_columns
        sushi_app_entry.next_dataset_keys = sushi_app_instance.next_dataset.keys
        sushi_app_entry.save
      end
    end

    # delete check
    delete_apps = SushiApplication.all.map{|app| app.class_name} - sushi_apps.map{|app| app.gsub(/\.rb/,'').gsub(/\.sh/,'')}
    delete_apps.each do |class_name|
      app = SushiApplication.find_by_class_name(class_name)
      SushiApplication.delete(app)
    end

    @sushi_apps = SushiApplication.all.sort_by{|app| app.analysis_category}
    render action: "index"
  end
end
