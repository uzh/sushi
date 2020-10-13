namespace :lib do
  desc "Check ref_selector"
  task ref_selector: :environment do
    require 'sushi_fabric'
    require './lib/global_variables'
    include GlobalVariables

    Rails.cache = {}

    require 'pp'
    pp ref_selector
  end

end
