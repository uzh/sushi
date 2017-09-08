class AddRunNameOrderIdToDataSets < ActiveRecord::Migration
  def change
    add_column :data_sets, :run_name_order_id, :string
  end
end
