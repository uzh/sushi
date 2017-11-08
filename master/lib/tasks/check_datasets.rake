namespace :ds do

  desc "Check datasets parent ID"
  task star: :environment do
    puts ["project", "name", "#samples", "created_at"].join("\t")
    DataSet.all.each do |data_set|
      if data_set.name =~ /star/i
        puts [data_set.project.number, data_set.name, data_set.samples.length, data_set.created_at.to_s.split.first].join("\t")
      end
    end
  end
end
