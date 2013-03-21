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


class SampleController < ApplicationController
  def index
    @samples = Sample.all
  end

  def show
    bf = Bfabric.new(SushiFabric::Application.config.bfabric_user, 
                     SushiFabric::Application.config.bfabric_password)
    
    @sample_id = params[:id]
    @sample = bf.get_sample(@sample_id)
    
    @extracts = []
    e_ids = bf.get_extract_ids_for_sample(@sample_id.to_i)
    e_ids.each do |e_id|
      e = bf.get_extract(e_id)
      @extracts << e
    end
  end

  def edit
    @samples = Sample.all
    @sample = @samples[params[:id].to_i-1]
  end
end
