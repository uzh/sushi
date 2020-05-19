class CreateJobs < ActiveRecord::Migration[5.2]
  def change
    create_table :jobs do |t|
      t.integer :submit_job_id
      t.integer :input_dataset_id
      t.integer :next_dataset_id

      t.timestamps
    end
  end
end
