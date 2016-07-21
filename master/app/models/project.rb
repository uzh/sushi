class Project < ActiveRecord::Base
#  attr_accessible :number
  has_many :data_sets

  def saved?
    flag = false
    self.data_sets.each do |data_set|
      if data_set.saved?
        flag = true
        break
      end
    end
    flag
  end
end
