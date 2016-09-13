class AddBFabricIdToDataSets < ActiveRecord::Migration
  def change
    add_column :data_sets, :bfabric_id, :integer
  end
end
