class AddNfcoreAppsToDataSets < ActiveRecord::Migration[7.0]
  def change
    add_column :data_sets, :nfcore_apps, :text
  end
end
