require 'savon'
if SushiFabric::Application.config.fgcz?
require 'fgcz'
end
require 'csv' 
require 'sushi_fabric'
require 'SushiWrap'
require 'drb/drb'

# TODO: why also has to force it from here?
module Savon
  module SOAP
    class Response
      def to_xml
        http.body.force_encoding('UTF-8')
      end
    end
  end
end

class ApplicationController < ActionController::Base
  rescue_from DeviseLdapAuthenticatable::LdapException do |exception|
    render :text => exception, :status => 500
  end
  protect_from_forgery
  
  if SushiFabric::Application.config.fgcz?
    before_action :authenticate_user!
  end
  @@workflow_manager = DRbObject.new_with_uri(SushiFabric::WORKFLOW_MANAGER)

  def all_sushi_applications
    non_sushi_apps = ['SushiWrap.rb', 'optparse_ex.rb', 'global_variables.rb']
    lib_dir = File.expand_path('../../../lib', __FILE__)
    sushi_apps = Dir[File.join(lib_dir, '*.rb')].select{|script| !non_sushi_apps.include?(File.basename(script))}.to_a.map{|script| File.basename(script)}
    sushi_apps.concat Dir[File.join(lib_dir, '*.sh')].map{|script| File.basename(script)}
    sushi_apps.each do |script|
      if script =~ /\.rb/
        class_name = script.gsub(/\.rb/,'')
        require class_name
      elsif script =~ /\.sh/
        sushi_wrap = SushiWrap.new(File.join(lib_dir, script))
        sushi_wrap.define_class
      end
    end
    sushi_apps
  end
  def refresh_sushi_application

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
        sushi_app_entry.description = sushi_app_instance.description
        sushi_app_entry.save
      end
    end

    # delete check
    delete_apps = SushiApplication.all.map{|app| app.class_name} - sushi_apps.map{|app| app.gsub(/\.rb/,'').gsub(/\.sh/,'')}
    delete_apps.each do |class_name|
      app = SushiApplication.find_by_class_name(class_name)
      SushiApplication.delete(app)
    end
  end
  def runnable_application(data_set_headers, refresh = true)
    if refresh
      refresh_sushi_application
    end
    sushi_apps = SushiApplication.all.select do |app|
      (app.required_columns - data_set_headers.map{|colname| colname.to_s.gsub(/\[.+\]/,'').strip}).empty?
    end
  end
  def sample_path(data_set)
    paths = []
    data_set.samples.each do |sample|
      sample.to_hash.each do |header, file|
        if (header.tag?('File') or header.tag?('Link')) and file !~ /http/
          paths << File.dirname(file)
        end
      end
    end
    paths.uniq!
    if path = paths.first
      path.split('/')[0,2].join('/')
    end
  end
end
