class AddSushiAppNameToDataSets < ActiveRecord::Migration
  def change
    add_column :data_sets, :sushi_app_name, :string
  end
end
