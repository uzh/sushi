require 'savon'
if SushiFabric::Application.config.fgcz?
require 'fgcz'
end
require 'csv' 
require 'sushi_fabric'
require 'SushiWrap'
require 'drb/drb'
require 'tmpdir'

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
  # Temporarily disable LDAP exception handling
  # rescue_from DeviseLdapAuthenticatable::LdapException do |exception|
  #   render :plain => exception, :status => 500
  # end
  #  protect_from_forgery
  protect_from_forgery with: :exception

  after_action :flash_to_headers
  
  # Temporarily disable authentication for development
  # if SushiFabric::Application.config.fgcz?
  #   before_action :authenticate_user!
  # end
  @@sushi_server = eval(SushiFabric::Application.config.sushi_server_class).new

  def employee?
    view_context.employee?
  end
  def user_projects
    view_context.user_projects
  end
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
    NilClass.class_eval do
      def [](x)
        ''
      end
    end
    sushi_apps.select{|app| app =~ /\.rb$/}.each do |app|
      class_name = app.gsub(/\.rb/,'')
      updated_at = File.stat(File.join(lib_dir, app)).mtime
      unless sushi_app = SushiApplication.find_by_class_name(class_name) and updated_at < sushi_app.updated_at
        if sushi_app 
          load File.join(lib_dir, app)
        end
        begin
          sushi_app_instance = eval(class_name).new
          sushi_app_instance.instance_variable_set(:@dataset, {})
          sushi_app_instance.instance_variable_set(:@result_dir, '')
          sushi_app_entry = (sushi_app || SushiApplication.new)
          sushi_app_entry.class_name = class_name
          sushi_app_entry.analysis_category = sushi_app_instance.analysis_category
          sushi_app_entry.required_columns = sushi_app_instance.required_columns
          sushi_app_entry.next_dataset_keys = sushi_app_instance.next_dataset.keys
          sushi_app_entry.description = sushi_app_instance.description
          sushi_app_entry.employee = sushi_app_instance.employee
          sushi_app_entry.save
        rescue => err
          warn err
          warn "#{class_name} cannot be imported"
        end
      end
    end
    NilClass.class_eval do
      undef_method :[]
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
  def employee_apps
    apps = {}
    SushiApplication.all.select{|app| app.employee}.each do |app|
      apps[app.class_name.to_s] = true
    end
    apps
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
  def save_dataset_tsv_in_gstore(data_set, data_set_file_name="dataset.tsv")
    if SushiFabric::Application.config.fgcz?
      target_dataset_tsv = ''
      Dir.mktmpdir do |dir|
        out_tsv = File.join(dir, data_set_file_name)
        data_set.save_as_tsv(out_tsv)
        project_number = session[:project]
        project = "p#{project_number}"
        dataset_path = if data_set.child
                         File.join(project, data_set.name)
                       elsif dirs = data_set.paths
                         if dirs.length > 1
                           File.join(project, data_set.name)
                         else
                           dirs.first
                         end
                       else
                         File.join(project, data_set.name)
                       end
        target_dir = File.join(SushiFabric::GSTORE_DIR, dataset_path)
        target_dataset_tsv = File.join(target_dir, data_set_file_name)
        #print File.read(out_tsv)
        commands = @@sushi_server.copy_commands(out_tsv, target_dir, "force")
        commands.each do |command|
          #puts command
          `#{command}`
        end
      end
    end
  end

  private

  def flash_to_headers
    return unless request.xhr?
    response.headers['X-Message'] = flash[:error] || flash[:alert] || flash[:notice]
    flash.discard
  end
end
