class AddFieldsToJobs < ActiveRecord::Migration[7.0]
  def change
    add_column :jobs, :stdout_path, :string
    add_column :jobs, :stderr_path, :string
    add_column :jobs, :submit_command, :text
    add_column :jobs, :status, :string
    add_column :jobs, :user, :string
    add_column :jobs, :start_time, :datetime
    add_column :jobs, :end_time, :datetime

    add_index :jobs, :status
  end
end
