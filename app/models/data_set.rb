class DataSet < ActiveRecord::Base
  attr_accessible :name, :md5
  has_many :samples
  belongs_to :project
  has_many :data_sets, :foreign_key => :parent_id
  belongs_to :data_set, :foreign_key => :parent_id

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
  def export_tsv(file_path)
    CSV.open(file_path, 'w', :col_sep=>"\t") do |out|
      out << headers
      self.samples.each do |sample|
        out << headers.map{|header| sample.to_hash[header]}
      end
    end
  end
end
