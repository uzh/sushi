class CreateDataSets < ActiveRecord::Migration[5.2]
  def change
    create_table :data_sets do |t|
      t.integer :project_id
      t.integer :parent_id
      t.string :name
      t.string :md5

      t.timestamps
    end
  end
end
