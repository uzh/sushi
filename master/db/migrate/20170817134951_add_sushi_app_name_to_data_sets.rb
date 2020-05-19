class AddSushiAppNameToDataSets < ActiveRecord::Migration[5.2]
  def change
    add_column :data_sets, :sushi_app_name, :string
  end
end
