require 'savon'
if `hostname`.chomp =~ /fgcz/
require 'fgcz'
end
require 'csv' 
require 'sushi_fabric'
require 'SushiWrap'

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
  
  if `hostname`.chomp =~ /fgcz/
    before_filter :authenticate_user!
  end

  def all_sushi_applications
    non_sushi_apps = ['SushiWrap.rb', 'optparse_ex.rb']
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
  def runnable_application(data_set_headers)
    sushi_apps = all_sushi_applications

    # filter application with data_set#required_columns
    lib_dir = File.expand_path('../../../lib', __FILE__)
    sushi_apps = sushi_apps.sort.select do |script|
      class_name = ''
      if script =~ /\.rb/
        class_name = script.gsub(/\.rb/,'')
        #require class_name
      elsif script =~ /\.sh/
        class_name = script.gsub(/\.sh/,'')
        #sushi_wrap = SushiWrap.new(File.join(lib_dir, script))
        #sushi_wrap.define_class
      end
      sushi_app = eval(class_name).new
      required_columns = sushi_app.required_columns
      (required_columns - data_set_headers.map{|colname| colname.to_s.gsub(/\[.+\]/,'').strip}).empty?
    end
    sushi_apps
  end
end
