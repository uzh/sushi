class AddRunnableAppsToDataSet < ActiveRecord::Migration
  def change
    add_column :data_sets, :runnable_apps, :text
  end
end
