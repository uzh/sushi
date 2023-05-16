class AddOrderIdToDataSets < ActiveRecord::Migration[6.1]
  def change
    add_column :data_sets, :order_id, :integer
    add_index :data_sets, :order_id
  end
end
