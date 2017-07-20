class AddDataSetTreeToProjects < ActiveRecord::Migration
  def change
    add_column :projects, :data_set_tree, :text
  end
end
