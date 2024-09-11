class Project < ActiveRecord::Base
#  attr_accessible :number
  has_many :data_sets
  serialize :data_set_tree, Hash

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
    register_command = "register_sushi_dataset_into_bfabric"
    check_command = "check_dataset_bfabric"
    if SushiFabric::Application.config.fgcz? and system("which #{register_command} > /dev/null 2>&1") and system("which #{check_command} > /dev/null 2>&1")
      t = Time.new(2016)
      self.data_sets.each do |data_set|
        if data_set.data_set.nil? and data_set.created_at >= t # if it is the top node dataset (== suppose raw dataset)
          data_set.register_bfabric(op)
        end
      end
    end
  end
  def make_tree_node(data_set)
      project_dataset_ids = Hash[*(self.data_sets.map{|data_set| [data_set.id, true]}.flatten)]
      node = {"id" => data_set.id,
              "text" => data_set.data_sets.length.to_s+" "+data_set.name+" <small><font color='gray'>"+data_set.comment.to_s+"</font></small>",
              "a_attr" => {"href"=>"/data_set/p#{self.number}/#{data_set.id}"}
              }
      if parent = data_set.data_set and project_dataset_ids[parent.id]
        node["parent"] = parent.id
      else
        node["parent"] = "#"
      end
      node
  end
  def construct_data_set_tree
    tree = {}
    project_dataset_ids = Hash[*(self.data_sets.map{|data_set| [data_set.id, true]}.flatten)]
    self.data_sets.each do |data_set|
      tree[data_set.id] = make_tree_node(data_set)
    end
    self.data_set_tree = tree
    self.save
  end
  def add_tree_node(data_set)
    node = make_tree_node(data_set)
    self.data_set_tree[data_set.id] = node
    self.save
  end
end
