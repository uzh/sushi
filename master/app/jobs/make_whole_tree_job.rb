class MakeWholeTreeJob < ApplicationJob
  queue_as :default

  def perform(project_id)
    # update whole tree in Redis
    workflow_manager = DRbObject.new_with_uri(SushiFabric::WORKFLOW_MANAGER)
    project = Project.find_by_id(project_id)
    root = []
    project_dataset_ids = Hash[*(project.data_sets.map{|data_set| [data_set.id, true]}.flatten)]
    project.data_sets.each do |data_set|
      node = {"id" => data_set.id,
              "text" => data_set.data_sets.length.to_s+" "+data_set.name+" <small><font color='gray'>"+data_set.comment.to_s+"</font></small>",
              "a_attr" => {"href"=>"/data_set/p#{project.number}/#{data_set.id}"}
              }
      if parent = data_set.data_set and project_dataset_ids[parent.id]
        node["parent"] = parent.id
      else
        node["parent"] = "#"
      end
      root << node
    end
    json = root.sort_by{|node| node["id"]}.reverse.to_json
    workflow_manager.save_dataset_tree(project.number, json)
  end
end
