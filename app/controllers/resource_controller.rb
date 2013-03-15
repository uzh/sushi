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

    new_samples = []
    if r and relativepath = r[:relativepath]
      sample = if hit = Sample.find_by_resource_id(params[:id])
                 hit
               else
                 new_sample = Sample.new
                 new_sample.path = "/srv/gstore/projects/"+r[:relativepath]
                 new_sample.name = if s
                                 s[:name]+" [e"+e[:@id]+"]"
                               else
                                 File.basename(sample.path)
                               end
                 new_sample.resource_id = params[:id].to_s
                 new_sample
               end
      new_samples << sample
    else
      result_paths = `public/wfm_get_result_paths #{params[:id].gsub(/^JobID/,'')}`
      result_paths = result_paths.split(/,/)
      result_paths.each_with_index do |path, i|
        search_id = if result_paths.length < 2
                      params[:id].to_s
                    else
                      params[:id].to_s + '_' + (i+1).to_s
                    end
        sample = if hit = Sample.find_by_resource_id(search_id)
                   hit
                 else
                   new_sample = Sample.new
                   new_sample.path = path
                   new_sample.name = File.basename(path)
                   new_sample.resource_id = if result_paths.length < 2
                                              params[:id].to_s
                                            else
                                              params[:id].to_s + '_' + (i+1).to_s
                                            end
                   new_sample
                 end
        new_samples << sample
      end
    end

    new_samples.each do |sample|
      session[:basket] << sample
    end
    
    redirect_to request.referer
  end
  
  def remove_from_basket
    if session[:basket] and params[:id]
      session[:basket].delete_if { |sample| sample.resource_id == params[:id].to_s }
    end
    redirect_to request.referer
  end
end
