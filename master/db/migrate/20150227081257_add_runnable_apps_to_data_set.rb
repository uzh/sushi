class AddRunnableAppsToDataSet < ActiveRecord::Migration[5.2]
  def change
    add_column :data_sets, :runnable_apps, :text
  end
end
