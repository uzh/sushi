class DataSet < ActiveRecord::Base
  attr_accessible :name, :parent_id
  has_many :samples
  has_many :data_sets, :foreign_key => :parent_id
  belongs_to :data_set, :foreign_key => :parent_id
end
