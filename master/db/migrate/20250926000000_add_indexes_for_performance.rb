class AddIndexesForPerformance < ActiveRecord::Migration[7.0]
  def up
    add_index :jobs, :next_dataset_id unless index_exists?(:jobs, :next_dataset_id)
    add_index :jobs, :input_dataset_id unless index_exists?(:jobs, :input_dataset_id)
    add_index :jobs, [:next_dataset_id, :id] unless index_exists?(:jobs, [:next_dataset_id, :id])

    add_index :data_sets, :project_id unless index_exists?(:data_sets, :project_id)
    add_index :data_sets, :parent_id unless index_exists?(:data_sets, :parent_id)

    add_index :projects, :number unless index_exists?(:projects, :number)

    add_index :samples, :data_set_id unless index_exists?(:samples, :data_set_id)
  end

  def down
    remove_index :jobs, column: [:next_dataset_id, :id] if index_exists?(:jobs, [:next_dataset_id, :id])
    remove_index :jobs, column: :input_dataset_id if index_exists?(:jobs, :input_dataset_id)
    remove_index :jobs, column: :next_dataset_id if index_exists?(:jobs, :next_dataset_id)

    remove_index :data_sets, column: :project_id if index_exists?(:data_sets, :project_id)
    remove_index :data_sets, column: :parent_id if index_exists?(:data_sets, :parent_id)

    remove_index :projects, column: :number if index_exists?(:projects, :number)

    remove_index :samples, column: :data_set_id if index_exists?(:samples, :data_set_id)
  end
end


