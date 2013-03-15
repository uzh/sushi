class Sample < ActiveRecord::Base
  attr_accessible :name, :path, :parent_id
  has_many :data_lists
  has_many :samples, :foreign_key => "parent_id"
  belongs_to :sample, :foreign_key => "parent_id"
end
