class AddWorkunitIdToDataSets < ActiveRecord::Migration[5.2]
  def change
    add_column :data_sets, :workunit_id, :integer
  end
end
