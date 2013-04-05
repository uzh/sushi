class CreateSamples < ActiveRecord::Migration
  def change
    create_table :samples do |t|
      t.integer :data_set_id
      t.string :key_value

      t.timestamps
    end
  end
end
