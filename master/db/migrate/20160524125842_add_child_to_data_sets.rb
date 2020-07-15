class AddChildToDataSets < ActiveRecord::Migration[5.2]
  def change
    add_column :data_sets, :child, :boolean, default: false, null: false
  end
end
