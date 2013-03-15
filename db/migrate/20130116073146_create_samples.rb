class CreateSamples < ActiveRecord::Migration
  def change
    create_table :samples do |t|
      t.string :name
      t.string :path
      t.integer :parent_id

      t.timestamps
    end
  end
end
