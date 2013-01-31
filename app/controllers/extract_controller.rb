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


class ExtractController < ApplicationController
  def show
    bf = Bfabric.new(SushiFabric::Application.config.bfabric_user, 
                     SushiFabric::Application.config.bfabric_password)
    
    @extract_id = params[:id]
    @extract = bf.get_extract(@extract_id)
    
    @resources = []
    r_ids = bf.get_resource_ids_for_extract(@extract_id.to_i)
    r_ids.each do |r_id|
      r = bf.get_resource(r_id)
      @resources << r
    end
  end
end
