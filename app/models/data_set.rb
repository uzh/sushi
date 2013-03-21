class DataSet < ActiveRecord::Base
  attr_accessible :name, :note, :parent_id
  has_many :data_lists
  has_many :data_sets, :foreign_key => :parent_id
  belongs_to :data_set, :foreign_key => :parent_id
end
