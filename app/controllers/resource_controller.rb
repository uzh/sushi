class ResourceController < ApplicationController
  def add_to_basket
    if ! session.include? :basket
      session[:basket] = Array.new()
    end
    if ! session[:basket].include? params[:id]
      session[:basket] << params[:id]
    end
    
    redirect_to request.referer
  end
end
