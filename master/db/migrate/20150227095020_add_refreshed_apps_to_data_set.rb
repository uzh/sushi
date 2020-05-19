class AddRefreshedAppsToDataSet < ActiveRecord::Migration[5.2]
  def change
    add_column :data_sets, :refreshed_apps, :boolean
  end
end
