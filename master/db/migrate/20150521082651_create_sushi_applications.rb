class CreateSushiApplications < ActiveRecord::Migration[5.2]
  def change
    create_table :sushi_applications do |t|
      t.string :class_name
      t.string :analysis_category
      t.text :required_columns
      t.text :next_dataset_keys

      t.timestamps
    end
  end
end
