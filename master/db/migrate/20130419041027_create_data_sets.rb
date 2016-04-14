class CreateDataSets < ActiveRecord::Migration
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
