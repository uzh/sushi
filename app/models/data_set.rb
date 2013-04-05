class DataSet < ActiveRecord::Base
  attr_accessible :name, :parent_id
  has_many :samples
  has_many :data_sets, :foreign_key => :parent_id
  belongs_to :data_set, :foreign_key => :parent_id

  def headers
    self.samples.map{|sample| sample.to_hash.keys}.flatten.uniq
  end
  def saved?
    flag = false
    self.samples.each do |sample|
      if sample.saved?
        flag = true
        break
      end
    end
    flag
  end
  def to_csv
    self.samples.each do |sample|
      headers.each do |header|
      end
    end
  end

end
