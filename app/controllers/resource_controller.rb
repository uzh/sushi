class ResourceController < ApplicationController
  def add_to_basket
    bf = Bfabric.new(SushiFabric::Application.config.bfabric_user,
                     SushiFabric::Application.config.bfabric_password)
    
    session[:basket] ||= []

    if session[:basket].include? params[:id]
      redirect_to request.referer
      return
    end


    r = bf.get_resource params[:id].to_i
    e = bf.get_extract r[:extract][:@id].to_i
    s = bf.get_sample e[:sample][:@id].to_i
    
    sample = if hit = Sample.find_by_resource_id(params[:id])
               hit
             else
               new_sample = Sample.new
               new_sample.name = s[:name]+" [e"+e[:@id]+"]"
               new_sample.path = "/srv/gstore/projects/"+r[:relativepath]
               new_sample.resource_id = params[:id].to_i
               new_sample
             end
    
    session[:basket] << sample
    
    redirect_to request.referer
  end
  
  def remove_from_basket
    if session[:basket] and params[:id]
      session[:basket].delete_if { |sample| sample.resource_id == params[:id].to_i }
    end
    redirect_to request.referer
  end
end
