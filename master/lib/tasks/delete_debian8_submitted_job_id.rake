namespace :job do

  desc "delete_debian8_submitted_job_id"
  task delete_debian8_submitted_job_id: :environment do
    File.readlines("lib/tasks/debian8_job_list_20200813.txt").to_a.map{|line| line.chomp}.each do |job_id|
      if job = Job.find_by_id(job_id)
        p job_id
        job.submitted_job_id = nil
        job.save
      end

    end
  end
end

