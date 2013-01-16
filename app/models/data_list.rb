class DataList < ActiveRecord::Base
  attr_accessible :data_set_id, :sample_id
  belongs_to :data_set
  belongs_to :sample
end
