class AddOrderIDsToDataSet < ActiveRecord::Migration[5.2]
  def change
    add_column :data_sets, :order_ids, :text
  end
end
