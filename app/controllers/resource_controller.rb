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
    e = if r and extract = r[:extract]
          bf.get_extract extract[:@id].to_i
        end
    s = if e and sample = e[:sample]
          bf.get_sample sample[:@id].to_i
        end
    
    sample = if hit = Sample.find_by_resource_id(params[:id])
               hit
             else
               new_sample = Sample.new
               new_sample.path = if r and relativepath = r[:relativepath]
                                   "/srv/gstore/projects/"+r[:relativepath]
                                 else
                                   `public/wfm_get_result_paths #{params[:id].gsub(/^9999/,'')}`
                                 end
               new_sample.name = if s
                                   s[:name]+" [e"+e[:@id]+"]"
                                 else
                                   File.basename(new_sample.path)
                                 end
               new_sample.resource_id = params[:id].to_s
               new_sample
             end
    
    session[:basket] << sample
    
    redirect_to request.referer
  end
  
  def remove_from_basket
    if session[:basket] and params[:id]
      session[:basket].delete_if { |sample| sample.resource_id == params[:id].to_s }
    end
    redirect_to request.referer
  end
end
