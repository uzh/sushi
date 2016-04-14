class AddUserIdToDataSet < ActiveRecord::Migration
  def change
    add_column :data_sets, :user_id, :integer
  end
end
