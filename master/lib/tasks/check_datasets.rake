namespace :ds do

  desc "Check datasets"
  task star: :environment do
    puts ["project", "name", "#samples", "created_at"].join("\t")
    DataSet.all.each do |data_set|
      if data_set.name =~ /star/i
        puts [data_set.project.number, data_set.name, data_set.samples.length, data_set.created_at.to_s.split.first].join("\t")
      end
    end
  end

  task roots: :environment do
    ds2nc = {}
    DataSet.all.each do |data_set|
      unless data_set.data_set
        # root data_set
        child_count = 0
        data_set.data_sets.each do |child|
          if child.name =~ /Fastqc/ or child.name =~ /FastqScreen/
            # nothing
          else
            child_count += 1
          end
        end
        ds2nc[data_set] = child_count
      end
    end
    puts ["project", "name", "#samples", "created_at", "link"].join("\t")
    ds2nc.sort_by{|data_set, child_count| data_set.project.number}.each do |data_set, child_count|
      link = "https://fgcz-sushi.uzh.ch/data_set/#{data_set.id}"
      puts [data_set.project.number, data_set.name, child_count, data_set.created_at.to_s.split.first, link].join("\t")
    end
  end

end
