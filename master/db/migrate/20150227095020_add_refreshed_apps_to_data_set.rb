class AddRefreshedAppsToDataSet < ActiveRecord::Migration
  def change
    add_column :data_sets, :refreshed_apps, :boolean
  end
end
