class Job < ActiveRecord::Base
#  attr_accessible :submit_job_id
  belongs_to :data_set, :foreign_key => :next_dataset_id
end
