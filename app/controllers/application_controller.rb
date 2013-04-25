require 'savon'
require 'bfabric'
require 'csv'
require 'sushiApp'

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
  
  before_filter :authenticate_user!  
  
  def get_user_id
    if !session[:bf_user_id]
      if !current_user
        logger.error "No logged in user."
        0
      else
        bf = Bfabric.new(SushiFabric::Application.config.bfabric_user,
                         SushiFabric::Application.config.bfabric_password)
        session[:bf_user_id] = bf.get_user_id current_user.login
      end
    end
    session[:bf_user_id]
  end
end
