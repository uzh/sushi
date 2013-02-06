class JobMonitoringController < ApplicationController
  def index
    @job_list = `public/wfm_job_list`
    @job_list = @job_list.split(/\n/).map{|job| job.split(/,/)}
  end
end
