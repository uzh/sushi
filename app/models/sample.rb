class Sample < ActiveRecord::Base
  attr_accessible :name, :path
  has_many :data_list
end
