class AddChildToDataSets < ActiveRecord::Migration
  def change
    add_column :data_sets, :child, :boolean, default: false, null: false
  end
end
