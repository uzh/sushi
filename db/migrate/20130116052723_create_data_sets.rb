class CreateDataSets < ActiveRecord::Migration
  def change
    create_table :data_sets do |t|
      t.string :note
      t.string :name
      t.integer :parent_id

      t.timestamps
    end
  end
end
