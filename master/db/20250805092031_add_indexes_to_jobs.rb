class AddIndexesToJobs < ActiveRecord::Migration[7.0]
  def change
    # Index for JOIN operations with next_dataset_id
    add_index :jobs, :next_dataset_id
    
    # Composite index for ORDER BY id DESC with JOIN optimization
    add_index :jobs, [:id, :next_dataset_id], order: { id: :desc }
    
    # Index for potential future use with input_dataset_id
    add_index :jobs, :input_dataset_id
  end
end
