class SampleController < ApplicationController
  def index
    @samples = Sample.all
  end
end
