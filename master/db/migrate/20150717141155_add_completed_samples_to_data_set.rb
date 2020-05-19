class AddCompletedSamplesToDataSet < ActiveRecord::Migration[5.2]
  def change
    add_column :data_sets, :completed_samples, :integer
  end
end
