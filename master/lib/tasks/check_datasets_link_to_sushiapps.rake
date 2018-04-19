namespace :ds do
  desc "Check all datases with a sushiapp as a argument"
  task show_all: :environment do
    sushi_app_name = ENV['SUSHI_APP_NAME']
    p DataSet.all.first
    p sushi_app_name
    count = 0
    DataSet.all.select{|data_set| 
      data_set.sushi_app_name =~ /#{sushi_app_name}/i
    }.sort_by{|data_set| 
      data_set.created_at
    }.reverse.each do |data_set|
      count += 1
      p data_set.sushi_app_name
      p data_set.created_at
    end
    p count
  end
end
