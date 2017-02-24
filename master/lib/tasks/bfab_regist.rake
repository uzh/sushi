namespace :bfab_regist do

  desc "Registration DataSet to B-Fabric"
  task run: :environment do
    project_number = ENV['project_number']
    puts "project_number: #{project_number}"
    if project = Project.find_by_number(project_number.to_i)
      project.register_bfabric
    else
      puts "There is no such project: p#{project_number.to_i}"
    end
  end

  task check_importable_datasets: :environment do
    total_registered_data_sets = 0
    t = Time.new(2016)
    Project.all.sort_by{|project| project.number}.each do |project|
      data_set_count = 0
      project.data_sets.each do |data_set| 
        if data_set.created_at >= t
          data_set_count += 1
        end
      end
      total_registered_data_sets += data_set_count
      puts "project_number: #{project.number}, data_sets: #{data_set_count}"
      #if project.number == 1535
      #  project.register_bfabric
      #end
      project.register_bfabric
    end
    puts "total_registered_data_sets: #{total_registered_data_sets}"
  end


  task count_importable_datasets: :environment do
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
