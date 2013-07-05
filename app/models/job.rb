class Job < ActiveRecord::Base
  attr_accessible :submit_job_id
  belongs_to :next_dataset, :class_name => "DataSet", :foreign_key => :next_dataset_id
end
