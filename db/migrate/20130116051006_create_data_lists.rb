class CreateDataLists < ActiveRecord::Migration
  def change
    create_table :data_lists do |t|
      t.integer :data_set_id
      t.integer :sample_id

      t.timestamps
    end
  end
end
