namespace :job do

  desc "Get all Job id in Debian8"
  task get_all: :environment  do
    Job.all.each do |job|
      puts job.id
    end
  end

end
