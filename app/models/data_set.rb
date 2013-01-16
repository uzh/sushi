class DataSet < ActiveRecord::Base
  attr_accessible :note
  has_many :data_lists
end
