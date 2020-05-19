class AddBFabricIdToDataSets < ActiveRecord::Migration[5.2]
  def change
    add_column :data_sets, :bfabric_id, :integer
  end
end
