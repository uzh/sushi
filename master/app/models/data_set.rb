class DataSet < ActiveRecord::Base
#  attr_accessible :name, :md5, :comment
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
  def factor_first_headers
    headers.sort_by do |col| 
      if col == 'Name' 
        0
      elsif col.scan(/\[(.*)\]/).flatten.join =~ /Factor/
        1
      else
        2
      end
    end
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
        out << headers.map{|header| 
          val = sample.to_hash[header]
          val.to_s.empty? ? nil:val}
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
  def register_bfabric(op = 'new')
    base = "public/register_sushi_dataset_into_bfabric"
    check = "public/check_dataset_bfabric"
    if SushiFabric::Application.config.fgcz? and File.exist?(base) and File.exist?(check)
      dataset_tsv = File.join(SushiFabric::Application.config.scratch_dir, "dataset.tsv")
      open(dataset_tsv, "w") do |out|
        out.print self.tsv_string
      end 
      puts "# created: #{dataset_tsv}"
      option_check = if op == 'new' and !self.bfabric_id
                       true
                     elsif op == 'update' and bfabric_id = self.bfabric_id
                       com = "#{check} #{bfabric_id}"
                       if out = `#{com}`
                         eval(out.chomp.downcase)
                       end
                     elsif op == 'renewal'
                       true
                     end
      command = if dataset_tsv and File.exist?(dataset_tsv)
                  if parent_dataset = self.data_set
                    if bfabric_id = parent_dataset.bfabric_id
                      [base, "p#{self.project.number}", dataset_tsv, self.name, self.id, bfabric_id].join(" ")
                    end
                  else
                    [base, "p#{self.project.number}", dataset_tsv, self.name, self.id].join(" ")
                  end
                end
      if option_check and command and bfabric_id = `#{command}`
        puts "$ #{command}"
        puts "# mode: #{op}"
        puts "bfabric_id: #{bfabric_id}"
        if bfabric_id.split(/\n/).uniq.length < 2 and bfabric_id.chomp.to_i > 0
          self.bfabric_id = bfabric_id.chomp.to_i
          puts "# BFabricID: #{self.bfabric_id}"
          self.save
        else
          puts "# Not executed properly:"
          puts "# BFabricID: #{bfabric_id}"
        end
        File.unlink dataset_tsv
        puts "# removed: #{dataset_tsv}"
      end
      if child_data_sets = self.data_sets
        child_data_sets.each do |child_data_set|
          child_data_set.register_bfabric(op)
        end
      end
    end
  end
end
