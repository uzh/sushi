namespace :data_delivery do

  desc "Check datasets only data delivery"
  task check: :environment do
    delivery_data_candidates = []
    DataSet.all.each do |data_set|
      unless data_set.data_set
        # root data_set
        flag = false
        data_set.data_sets.each do |child|
          if child.name =~ /Fastqc/ or child.name =~ /FastqScreen/
            # nothing
          else
            flag = true
            break
          end
        end
        unless flag
          delivery_data_candidates << data_set
        end
      end
    end
    puts ["project", "name", "#children", "created_at"].join("\t")
    delivery_data_candidates.each do |data_set|
      puts [data_set.project.number, data_set.name, data_set.data_sets.length, data_set.created_at.to_s.split.first].join("\t")
    end
    # check
    # p delivery_data_candidates.length
    # 980
  end
end
