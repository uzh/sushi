class CreateSamples < ActiveRecord::Migration
  def change
    create_table :samples do |t|
      t.string :key_value

      t.timestamps
    end
  end
end
