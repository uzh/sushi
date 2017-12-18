namespace :update_DB_ExtractID2NewSampleID do

  desc "Replace current ExtractID to new SampleID for new release of BFabric"
  task run: :environment do
    require 'csv'
    all_samples_csv = ENV['ALL_SAMPLES_CSV']
    warn all_samples_csv
    eid2newid = {}
    sid2newid = {}
    edup_cases = {}
    sdup_cases = {}
    CSV.table(all_samples_csv, headers: true).each do |row|
      eid = row[:oldextractid]
      newid = row[:id]
      unless eid2newid[eid]
        eid2newid[eid] = newid
      else
        unless eid.to_i == 0
          edup_cases[eid] ||= []
          edup_cases[eid] << newid
        end
      end

      sid = row[:oldsampleid]
      unless sid2newid[sid]
        sid2newid[sid] = newid
      else
        unless sid.to_i == 0
          sdup_cases[sid] ||= []
          sdup_cases[sid] << newid
        end
      end

    end
    edup_cases.each do |eid, newids|
      warn "Old ExtractID: #{eid}"
      newids.each do |id|
        warn "\tNew SampleID: #{id}"
      end
    end
    sdup_cases.each do |sid, newids|
      warn "Old SampleID: #{sid}"
      newids.each do |id|
        warn "\tNew SampleID: #{id}"
      end
    end

    warn
    warn "# old extract dup case: #{edup_cases.keys.length}"
    warn "# old extract total: #{eid2newid.keys.length}"
    warn "# old sample  dup case: #{sdup_cases.keys.length}"
    warn "# old sample  total: #{sid2newid.keys.length}"
    warn

    total_datasets = DataSet.all.length
    warn "# total datasets: #{total_datasets}"
    DataSet.all.each.with_index do |data_set, i|
      if i%100==0
        warn "#{i}/#{total_datasets}"
      end
      headers = data_set.headers
      if headers.any?{|header| header =~ /Extract Id/}
        data_set.samples.each do |sample|
          sample_hash = sample.to_hash
          if extract_id_ = (sample_hash["Extract Id [B-Fabric]"]||sample_hash["Extract Id"]) and !extract_id_.to_s.empty? and extract_id_.to_s != "NA" and extract_id_.to_s != "unk" and extract_id_.to_s != "?" 
            extract_id = extract_id_.gsub(/bfe_/, '')
            if newid = eid2newid[extract_id.to_i]
              puts sample_hash
              puts "data_set.id: #{data_set.id}"
              puts "sample.id: #{sample.id}"
              puts "Old Extract Id: #{extract_id}"

              sample_hash["Sample Id [B-Fabric]"] = "bfs_#{newid}"
              sample_hash.delete("Extract Id [B-Fabric]")
              sample.key_value = sample_hash.to_s
              sample.save
            else
              warn sample_hash
              warn "data_set.id: #{data_set.id}"
              warn "sample.id: #{sample.id}"
              warn "Old Extract Id: #{extract_id}"
              warn "No corresponding new sample id"

              sample_hash.delete("Extract Id [B-Fabric]")
              sample_hash.delete("Extract Id")
              sample.key_value = sample_hash.to_s
              sample.save
            end
          else
              warn sample_hash
              warn "data_set.id: #{data_set.id}"
              warn "sample.id: #{sample.id}"
              warn "Old Extract Id: #{extract_id}"
              warn "No corresponding new sample id"


              sample_hash.delete("Extract Id [B-Fabric]")
              sample_hash.delete("Extract Id")
              sample.key_value = sample_hash.to_s
              sample.save
          end
        end
      end
    end
  end
end
