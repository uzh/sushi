class Sample < ActiveRecord::Base
  attr_accessible :name, :path, :parent_id
  has_many :data_lists
end
