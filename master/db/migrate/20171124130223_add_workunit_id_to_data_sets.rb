class AddWorkunitIdToDataSets < ActiveRecord::Migration
  def change
    add_column :data_sets, :workunit_id, :integer
  end
end
