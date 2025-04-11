namespace :ds do

  def sample_path(data_set)
    paths = []
    data_set.samples.each do |sample|
      sample.to_hash.each do |header, file|
        if (header.tag?('File') or header.tag?('Link')) and file !~ /http/
          paths << File.dirname(file)
        end
      end
    end
    paths.uniq!
    paths
  end

  desc "Check STAR BWA Bowtie datasets"
  task star_bwa_bowtie: :environment do
    puts ["project", "name", "#samples", "created_at", "link", "gstore_path", "#bams", "size [GB]"].join("\t")
    DataSet.all.each do |data_set|
      if data_set.name =~ /star/i or data_set.name =~ /bwa/i or data_set.name =~ /bowtie/i
        link = "https://fgcz-sushi.uzh.ch/data_set/p#{data_set.project.number}/#{data_set.id}"
        paths = sample_path(data_set)
        paths.delete('.')
        paths = paths.uniq.compact
        bams_size_total = 0
        dir_size_total = 0
        gstore_paths = []
        paths.each do |path|
          gstore_path = File.join("/srv/gstore/projects", path)
          gstore_paths << gstore_path
          bams_size = Dir[File.join(gstore_path, "*.bam")].to_a.length
          bams_size_total += bams_size
          com = "du -s #{gstore_path}"
          dir_size = if File.exist?(gstore_path)
                       `#{com}`.to_i
                     else
                       0
                     end
          dir_size_total += dir_size
        end
        dir_size_total = "%d" % (dir_size_total/1000000.0)
        puts [data_set.project.number, data_set.name, data_set.samples.length, data_set.created_at.to_s.split.first, link, gstore_paths.join(","), bams_size_total, dir_size_total].join("\t")
      end
    end
  end

  desc "Check root datasets"
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

  desc "Check only data delivery datasets"
  task data_delivery_candidates: :environment do
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

  desc "Save all datasets with date and SUSHIApp"
  task all_dataset_with_date_and_sushiapp: :environment do
    out_file = "all_datasets_with_sushiapp_#{Time.now.strftime("%Y-%m-%d")}.csv"
    headers = ["name", "project", "date", "SUSHIApp"]
    require 'csv'
    sushiapp2date = {}
    csv = CSV.generate("", headers: headers, write_headers: true, col_sep: ",") do |out|
      DataSet.all.each do |dataset|
        out << [dataset.name, dataset.project.number, dataset.updated_at.strftime("%Y-%m-%d"), dataset.sushi_app_name]
        if dataset.sushi_app_name
          sushiapp2date[dataset.sushi_app_name] ||= []
          sushiapp2date[dataset.sushi_app_name] << dataset.updated_at
        end
      end
    end
    File.write(out_file, csv)
    warn "# #{out_file} generated"

    out_file = "sushiapp_call_times.csv"
    headers = ["SUSHIApp", "call times", "latest"]
    csv = CSV.generate("", headers: headers, write_headers: true, col_sep: ",") do |out|
      sushiapp2date.keys.sort.each do |sushiapp|
        out << [sushiapp, sushiapp2date[sushiapp].length, sushiapp2date[sushiapp].max.strftime("%Y-%m-%d")]
      end
    end
    File.write(out_file, csv)
    warn "# #{out_file} generated"
  end

  desc "Save all datasets (samples)1"
  task dump_all_datasets_samples: :environment do
    out_file = "all_datasets_smaples_#{Time.now.strftime("%Y-%m-%d")}.csv"
    headers = ["name", "project", "samples"]
    require 'csv'
    csv = CSV.generate("", headers: headers, write_headers: true, col_sep: ",") do |out|
      DataSet.all.each do |dataset|
        row = [dataset.name, dataset.project.number]
        row << dataset.samples.map{|sample| sample.key_value}.join(":")
        out << row
      end
    end
    File.write(out_file, csv)
    warn "# #{out_file} generated"
  end

  desc "Save all datasets (samples)2"
  task dump_all_datasets_samples2: :environment do
    out_file = "all_datasets_smaples2_#{Time.now.strftime("%Y-%m-%d")}.txt"
    res = ActiveRecord::Base.connection.select_all("select * from samples;")
    File.write(out_file, res.to_a.join(":"))
    warn "# #{out_file} generated"
  end

  desc "Initialize DataSet.order_id"
  task init_order_ids: :environment do
    DataSet.order("id").each_with_index do |dataset, i|
      if dataset.data_set.nil? and dataset.order_ids.uniq.length == 1
        order_id = dataset.order_ids.first.to_i
        dataset.order_id = order_id
        dataset.save
        p [i, dataset.id, dataset.order_ids, dataset.order_id].join(",")
      end
    end
  end

  desc "All dataset and read counts"
  task all_dataset_read_counts: :environment do
    headers = ["DataSetID", "Name", "ReadCounts", "Date"]
    puts headers.join("\t")
    DataSet.order("id").each_with_index do |dataset, i|
      #if dataset.data_set.nil? and dataset.order_ids.uniq.length == 1
      if dataset.data_set.nil?
        read_counts = 0
        dataset.samples.each do |sample|
          if read_count = sample.to_hash["Read Count"]
            read_counts += read_count.to_i
          end
        end
        date = dataset.created_at.strftime("%Y-%m-%d")
        puts [dataset.id, dataset.name, read_counts, date].join("\t")
      end
    end
  end

  desc "Check all datasets with registered in 3 days"
  task check_unregistered_datasets: :environment do
    t0 = Time.now
    today = Date.today
    #p today
    #p (today - 3)
    #p today.to_s
    #p (today - 3).to_s
    unregistered_datasets = []
    unregistered_datasets_without_project = []
    #puts ["ID", "Date", "Project", "OrderIDs"].join("\t")
    puts ["ID", "Date", "Project"].join("\t")
    DataSet.order("id").each_with_index do |dataset, i|
        date = dataset.created_at.strftime("%Y-%m-%d")
        if date > (today-365*3).to_s
          unless dataset.bfabric_id
            order_ids = {}
            #dataset.samples.each do |sample|
            #  order_id = sample.to_hash["Order Id [B-Fabric]"]
            #  order_ids[order_id]
            #end
            #order_ids_ = unless order_ids.keys.empty?
            #               order_ids.keys.join(";")
            #             else
            #               ""
            #             end
            unregistered_datasets << dataset
            if project = dataset.project
              #puts [dataset.id, date.to_s, "p#{project.number}", order_ids_].join("\t")
              puts [dataset.id, date.to_s, "p#{project.number}"].join("\t")
            else
              unregistered_datasets_without_project << dataset
            end
          end
       end
    end
    puts "# #unregistered_datasets: #{unregistered_datasets.length}"
    puts "# #unregistered_datasets_without_project: #{unregistered_datasets_without_project.length}"
    puts "# run time: #{"%.2f" % (Time.now - t0)} [s]"
  end

  desc "Check all datasets in 2023"
  task check_all_datasets_in_2023: :environment do
    t0 = Time.now
    first_date = Date.new(2023,1,1)
    datasets = []
    puts ["ID", "Name", "Project", "SushiApp", "Samples", "Who", "Created", "BFabricID"].join("\t")
    DataSet.order("id").each_with_index do |dataset, i|
        date = dataset.created_at
        user = if user = dataset.user
                 user.login
               else
                 "sushi_lover"
               end
        if date > first_date and dataset.project
          datasets << dataset
          puts [dataset.id, dataset.name, dataset.project.number, dataset.sushi_app_name.to_s, "#{dataset.completed_samples.to_i}/#{dataset.num_samples.to_i}", user, date.strftime("%Y-%m-%d"), dataset.bfabric_id.to_s].join("\t")
       end
    end
    warn "# #datasets: #{datasets.length}"
    warn "# run time: #{"%.2f" % (Time.now - t0)} [s]"
  end

  desc "Count #sushiapps in datasets and make sushiapp ranks"
  task check_sushiapps_in_datasets_2023: :environment do
    t0 = Time.now
    first_date = Date.new(2023,1,1)
    sushiapp2count = {}
    DataSet.order("id").each_with_index do |dataset, i|
        date = dataset.created_at
        user = if user = dataset.user
                 user.login
               else
                 "sushi_lover"
               end
        if date > first_date and dataset.project and sushiapp = dataset.sushi_app_name
          sushiapp2count[sushiapp] ||= 0
          sushiapp2count[sushiapp] += 1
          #puts [dataset.id, dataset.name, dataset.project.number, dataset.sushi_app_name.to_s, "#{dataset.completed_samples.to_i}/#{dataset.num_samples.to_i}", user, date.strftime("%Y-%m-%d"), dataset.bfabric_id.to_s].join("\t")
       end
    end
    noused_sushiapps = SushiApplication.all.to_a.map{|sushiapp| sushiapp.class_name}
    puts ["SUSHIApp", "#Used"].join("\t")
    sushiapp2count.sort_by{|sushi_app_name, count| count}.reverse.each do |sushi_app_name, count|
      puts [sushi_app_name, count].join("\t")
      noused_sushiapps.delete(sushi_app_name)
    end
    noused_sushiapps.sort.each do |sushi_app_name|
      puts [sushi_app_name, "0"].join("\t")
    end
    warn "# #sushiapps: #{sushiapp2count.keys.length}"
    warn "# run time: #{"%.2f" % (Time.now - t0)} [s]"
  end

  def dfs_sort_datasets(root_dataset, selected_datasets, result = [])
    if selected_datasets.include?(root_dataset)
      result << root_dataset
    end

    root_dataset.data_sets.order(:id).each do |child|
      dfs_sort_datasets(child, selected_datasets, result)
    end

    result
  end
  def out_dataset_list(out_file, datasets)
    class << datasets
      attr_accessor :samples
    end
    samples = datasets.inject(0){|sum, dataset| sum + dataset.completed_samples.to_i}
    datasets.samples = samples
    open(out_file, "w") do |out|
      out.puts ["ID", "Project", "Name", "SushiApp", "Samples", "Who", "Created", "BFabricID", "Order IDs"].join("\t") 
      datasets.each do |dataset|
        date = dataset.created_at
        user = if user = dataset.user
                 user.login
               else
                 "sushi_lover"
               end
        out.puts [dataset.id, dataset.project.number, dataset.name, dataset.sushi_app_name.to_s, "#{dataset.completed_samples.to_i}/#{dataset.num_samples.to_i}", user, date.strftime("%Y-%m-%d"), dataset.bfabric_id.to_s, dataset.order_ids.join(",")].join("\t")
      end
      out.puts "# #datasets: #{datasets.length}"
      out.puts "# #samples:  #{samples}"
    end
    warn "# [#{Time.now.strftime("%Y-%m-%d %H:%M:%S")}] #{out_file} generated"
  end

  def collect_all_child_datasets(parent_dataset, sorted_datasets, result = [])
    children = sorted_datasets.select { |dataset| dataset.parent_id == parent_dataset.id }
    result.concat(children)
    children.each do |child|
      collect_all_child_datasets(child, sorted_datasets, result)
    end
    result
  end
 
  desc "Register datasets to BFabric"
  task :register_datasets, [:year,:run] => :environment do |t, args|
    # bundle exec rake ds:register_datasets[2024] RAILS_ENV=production DISABLE_DATABASE_ENVIRONMENT_CHECK=1
    # bundle exec rake ds:register_datasets[2024,run] RAILS_ENV=production DISABLE_DATABASE_ENVIRONMENT_CHECK=1
    dataset_list_log = "selected_dataset_list.log"
    dataset_tree_log = "sorted_dataset_list.log"
    black_list_log = "black_dataset_list.log"

    run = args[:run]
    t0 = Time.now
    year = if year_ = args[:year]
             year_.to_i
           else
             Date.today.year
           end
    first_date = Date.new(year,1,1)
    selected_datasets = []

    # keep unsorted dataset list
    selected_datasets = DataSet.order("id").select do |dataset|
      date = dataset.created_at
      date > first_date && dataset.project && dataset.bfabric_id.nil?
    end
    out_dataset_list(dataset_list_log, selected_datasets)

    # check order_ids
    selected_datasets.each do |dataset|
      if dataset.order_ids.empty?
        dataset.check_order_ids
      end
    end

    # sorting dataset list
    sorted_datasets = []
    root_datasets = selected_datasets.select do |dataset|
      dataset.data_set.nil? || !selected_datasets.include?(dataset.data_set)
    end
    root_datasets.each do |root_dataset|
      dfs_sort_datasets(root_dataset, selected_datasets, sorted_datasets)
    end
    
    # keep child datasets
    child_datasets = {}
    sorted_datasets.each do |parent_dataset|
      child_datasets[parent_dataset] = collect_all_child_datasets(parent_dataset, sorted_datasets)
    end
    #child_datasets.each do |parent_dataset, children|
    #  puts "#{parent_dataset.id}: #{children.map{|child| child.id.to_s}.join(",")}"
    #end 

    black_list_datasets = []
    # filter out expected failure datasets and the children
    sorted_datasets.each do |dataset|
      if !black_list_datasets.include?(dataset) and (parent_dataset = dataset.data_set and (parent_dataset.bfabric_id.nil? or black_list_datasets.include?(parent_dataset))) or
        dataset.order_ids.empty? or
        dataset.completed_samples.to_i < dataset.num_samples.to_i or
        dataset.sushi_app_name.nil?

        black_list_datasets << dataset
        black_list_datasets.concat(child_datasets[dataset])
        #sorted_datasets -= [dataset].concat(child_datasets[dataset])
      end
    end
    black_list_datasets.uniq!
    sorted_datasets -= black_list_datasets
    out_dataset_list(black_list_log, black_list_datasets)
    out_dataset_list(dataset_tree_log, sorted_datasets)

    # for run
    if run
      puts ["dataset.id", "project.number", "dataset.name", "sushi_app_name", "completed_samples/num_samples", "user", "date", "bfabric_id", "order_ids"].join("\t")
      sorted_datasets.each do |dataset|
          date = dataset.created_at
          user = if user = dataset.user
                   user.login
                 else
                   "sushi_lover"
                 end
          puts [dataset.id, dataset.project.number, dataset.name, dataset.sushi_app_name.to_s, "#{dataset.completed_samples.to_i}/#{dataset.num_samples.to_i}", user, date.strftime("%Y-%m-%d"), dataset.bfabric_id.to_s, dataset.order_ids.join(",")].join("\t")
          dataset.register_bfabric
          sleep 1
          puts
      end
    end

    puts "# #selected_datasets: #{selected_datasets.length}"
    puts "# #selected_datasets.samples: #{selected_datasets.samples}"
    puts "# #sorted_datasets: #{sorted_datasets.length}"
    puts "# #sorted_datasets.samples: #{sorted_datasets.samples}"
    puts "# #black_list_datasets: #{black_list_datasets.length}"
    puts "# #black_list_datasets.samples: #{black_list_datasets.samples}"
    puts "# run time: #{"%.2f" % (Time.now - t0)} [s] (#{Time.now.strftime("%Y%m%d-%H:%M:%S")})"
  end
  class DataSet
    def search_order_id(recursive = 1)
      if parent = self.data_set 
        if !parent.order_ids.empty?
          print "\t"*recursive
          puts "-parent dataset ID:#{parent.id}-parent.order_ids:#{parent.order_ids.join(",")}"
          parent.order_ids.join(",")
        else
          print "\t"*recursive
          puts "-parent dataset ID:#{parent.id}"
          parent.search_order_id(recursive + 1)
        end
      end
    end
    def add_order_id_column(order_ids)
      # add new column
      new_header = "Order Id [B-Fabric]"
      self.samples.each_with_index do |sample, i|
        new_sample = sample.to_hash
        new_sample[new_header] = order_ids.strip
        self.samples[i].key_value = new_sample.to_s
        self.samples[i].save
      end
      self.md5 = self.md5hexdigest
      self.save
      #save_dataset_tsv_in_gstore(self) #=> error, not found
      puts "\t# Done: adding Order Ids (#{order_ids}) in dataset(#{self.id})"
    end
  end
  desc "No Order ID datasets list"
  task :no_order_id_datasets, [:year,:search_order_id,:update] => :environment do |t, args|
    # bundle exec rake ds:no_order_id_datasets[2024] RAILS_ENV=production DISABLE_DATABASE_ENVIRONMENT_CHECK=1
    search_order_id = args[:search_order_id]
    update = args[:update]

    t0 = Time.now
    year = if year_ = args[:year]
             year_.to_i
           else
             Date.today.year
           end
    first_date = Date.new(year,1,1)
    puts ["ID", "Name", "Project", "SushiApp", "Samples", "Who", "Created"].join("\t")
    datasets = []
    not_order_ids = 0
    add_order_ids = 0
    DataSet.order("id").each_with_index do |dataset, i|
      date = dataset.created_at
      user = if user = dataset.user
               user.login
             else
               "sushi_lover"
             end
      if date > first_date and dataset.project and dataset.bfabric_id.nil? and dataset.order_ids.empty?
        datasets << dataset
        puts [dataset.id, dataset.name, dataset.project.number, dataset.sushi_app_name.to_s, "#{dataset.completed_samples.to_i}/#{dataset.num_samples.to_i}", user, date.strftime("%Y-%m-%d")].join("\t")
        if search_order_id
          if order_ids = dataset.search_order_id
            puts "\t# Order ID: #{order_ids}"
            if update
              dataset.add_order_id_column(order_ids)
              add_order_ids += 1
            end
          else
            puts "\t# Order ID not found"
            not_order_ids += 1
          end
        end
      end
    end

    puts "# #datasets: #{datasets.length} (no Order Ids: #{not_order_ids}, added Order Ids: #{add_order_ids})"
    puts "# run time: #{"%.2f" % (Time.now - t0)} [s]"
  end

  desc "Check datasets with sushi_app_name missing"
  task check_datasets_sushi_app_name_missing: :environment do
    t0 = Time.now
    count_datasets_sushi_app_name_missing = 0
    count_updated = 0

    DataSet.all.each do |data_set|
      next if data_set.sushi_app_name

      count_datasets_sushi_app_name_missing += 1

      sample_path = data_set.sample_paths.first
      next unless sample_path

      parameters_tsv = File.join(SushiFabric::GSTORE_DIR, sample_path, "parameters.tsv")
      next unless File.exist?(parameters_tsv)

      # load parameters.tsv and search sushi_app_name
      sushi_app_name_value = nil
      File.foreach(parameters_tsv) do |line|
        key, value = line.chomp.split("\t", 2)
        if key == "sushi_app"
          sushi_app_name_value = value
          break
        end
      end

      # if sushi_app_name found, update data_set
      if sushi_app_name_value
        data_set.sushi_app_name = sushi_app_name_value
        data_set.save
        count_updated += 1
        puts "# Updated dataset #{data_set.id} with sushi_app_name: #{sushi_app_name_value}"
      end
    end

    puts "# datasets missing sushi_app_name: #{count_datasets_sushi_app_name_missing}"
    puts "# datasets updated: #{count_updated}"
    puts "# run time: #{"%.2f" % (Time.now - t0)} [s]"
  end

end

