class DataSet < ActiveRecord::Base
  attr_accessible :name
  has_many :samples

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
