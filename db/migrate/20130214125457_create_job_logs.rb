class CreateJobLogs < ActiveRecord::Migration
  def change
    create_table :job_logs do |t|
      t.integer :job_id
      t.string :result_link

      t.timestamps
    end
  end
end
