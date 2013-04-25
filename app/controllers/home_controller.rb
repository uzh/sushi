class HomeController < ApplicationController
  def index
    # tentatively only for develop
    session[:project] = 1001
  end
end
