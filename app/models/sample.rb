class Sample < ActiveRecord::Base
  attr_accessible :key_value
  belongs_to :data_set

  def to_hash
    eval(self.key_value)
  end
  def saved?
    if Sample.find_by_key_value(self.key_value)
      true
    else
      false
    end
  end
end
