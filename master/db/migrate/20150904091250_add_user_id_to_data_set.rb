class AddUserIdToDataSet < ActiveRecord::Migration[5.2]
  def change
    add_column :data_sets, :user_id, :integer
  end
end
