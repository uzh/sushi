class CreateSamples < ActiveRecord::Migration
  def change
    create_table :samples do |t|

      t.timestamps
    end
  end
end
