class HomeController < ApplicationController
  before_filter :authenticate_user!  
  def index
    if current_user
      redirect_to :user_root
      return
    end
  end
end
