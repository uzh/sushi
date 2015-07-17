class AddNumSamplesToDataSet < ActiveRecord::Migration
  def change
    add_column :data_sets, :num_samples, :integer
  end
end
