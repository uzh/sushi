require 'savon'
require 'bfabric'

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

class ProjectController < ApplicationController
  def index
    user_id = 1954
    bf = Bfabric.new(SushiFabric::Application.config.bfabric_user, 
                     SushiFabric::Application.config.bfabric_password)
    
    @projects = []
    p_ids = bf.get_project_ids_for_user(user_id) 
    p_ids.each do |p_id|
      p = bf.get_project(p_id)
      if p != nil
        @projects << p
      end
    end
  end
  
  def show
    bf = Bfabric.new(SushiFabric::Application.config.bfabric_user, 
                     SushiFabric::Application.config.bfabric_password)
    
    @project_id = params[:id]
    @samples = []
    s_ids = bf.get_sample_ids_for_project(@project_id.to_i)
    s_ids.each do |s_id|
      s = bf.get_sample(s_id)
      #s[:extracts] = []
      #e_ids = bf.get_extract_ids_for_sample(s_id.to_i)
      #e_ids.each do |e_id|
      #  e = bf.get_extract(e_id.to_i)
      #  s[:extracts] << e
      #end
      @samples << s
    end
  end
end
