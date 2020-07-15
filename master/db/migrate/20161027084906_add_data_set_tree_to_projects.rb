class AddDataSetTreeToProjects < ActiveRecord::Migration[5.2]
  def change
    add_column :projects, :data_set_tree, :text
  end
end
