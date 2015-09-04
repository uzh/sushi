class DataSet < ActiveRecord::Base
  attr_accessible :name, :md5, :comment
  has_many :samples
  has_many :jobs, :foreign_key => :next_dataset_id
  belongs_to :project
  has_many :data_sets, :foreign_key => :parent_id
  belongs_to :data_set, :foreign_key => :parent_id
  serialize :runnable_apps, Hash
  belongs_to :user

  def headers
    self.samples.map{|sample| sample.to_hash.keys}.flatten.uniq
  end
  def saved?
    if DataSet.find_by_md5(md5hexdigest)
      true
    else
      false
    end
  end
  def md5hexdigest
    key_value = self.samples.map{|sample| sample.key_value}.join + self.parent_id.to_s + self.project_id.to_s
    Digest::MD5.hexdigest(key_value)
  end
  def tsv_string
    string = CSV.generate(:col_sep=>"\t") do |out|
      out << headers
      self.samples.each do |sample|
        out << headers.map{|header| sample.to_hash[header]}
      end
    end
    string
  end
  def samples_length
    unless self.num_samples
      self.num_samples=self.samples.length
      self.save
    end
    self.num_samples
  end
end
