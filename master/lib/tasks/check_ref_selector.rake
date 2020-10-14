namespace :lib do
  desc "Check ref_selector"
  task ref_selector: :environment do
    require 'sushi_fabric'
    require './lib/global_variables'
    include GlobalVariables

    Rails.cache = {}

    require 'pp'
    pp ref_selector

    puts
    reference = "Eleusine_coracana/FGCZ/PR202_v2/Annotation/Release_01-2019-01-24"
    genome_root = fetch_genome_root(reference)
    p "reference=#{reference}"
    p "genome_root=#{genome_root}"
  end

end
