namespace :bfab_regist do

  desc "Registration DataSet to B-Fabric"
  task test: :environment do
    project_number = ENV['project_number']
    puts "project_number: #{project_number}"
    if project = Project.find_by_number(project_number.to_i)
      if data_sets = project.data_sets 
        t = Time.new(2016)
        data_sets.each do |data_set|
          if data_set.data_set.nil? and data_set.created_at >= t
            p [data_set.name, data_set.created_at]
          end
        end
      end
    else
      puts "There is no such project: p#{project_number.to_i}"
    end
  end

  task count_registable_datasets: :environment do
    t = Time.new(2016)
    counts = {}
    counts_total = 0
    DataSet.all.each do |data_set|
      if data_set.created_at >= t
        project_number = data_set.project.number
        counts[project_number] ||= 0
        counts[project_number] += 1
        counts_total += 1
      end
    end
    puts "total importable datasets : #{counts_total}"
  end
end
