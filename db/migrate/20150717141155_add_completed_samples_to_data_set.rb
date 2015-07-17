class AddCompletedSamplesToDataSet < ActiveRecord::Migration
  def change
    add_column :data_sets, :completed_samples, :integer
  end
end
