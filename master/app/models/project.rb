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
  def register_bfabric(op = 'new')
    base = "public/register_sushi_dataset_into_bfabric"
    check = "public/check_dataset_bfabric"
    if SushiFabric::Application.config.fgcz? and File.exist?(base) and File.exist?(check)
      t = Time.new(2016)
      self.data_sets.each do |data_set|
        if data_set.data_set.nil? and data_set.created_at >= t # if it is the top node dataset (== suppose raw dataset)
          data_set.register_bfabric(op)
        end
      end
    end
  end
end
