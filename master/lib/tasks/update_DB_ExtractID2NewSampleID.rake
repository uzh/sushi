namespace :update_DB_ExtractID2NewSampleID do

  desc "Replace current ExtractID to new SampleID for new release of BFabric"
  task run: :environment do
    require 'csv'
    all_samples_csv = ENV['ALL_SAMPLES_CSV']
    puts all_samples_csv
    eid2sid = {}
    dup_cases = {}
    CSV.foreach(all_samples_csv, headers: true) do |row|
      #puts "Old Extract Id: #{row["Old Extract Id"]}, New Sample Id: #{row["Sample Id"]}"
      eid = row["Old Extract Id"]
      sid = row["Sample Id"]

      unless eid2sid[eid]
        eid2sid[eid] = sid
      else
        unless eid.empty?
          dup_cases[eid] ||= []
          dup_cases[eid] << sid
        end
      end
    end
    dup_cases.each do |eid, sids|
      puts "ExtractID: #{eid}"
      sids.each do |sid|
        puts "\tNew SampleID: #{sid}"
      end
    end
    puts
    puts "# dup case: #{dup_cases.keys.length}"
    puts "# total: #{eid2sid.keys.length}"
  end
end
