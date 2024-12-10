class AddUserAndStartTimeAndEndTimeToJobs < ActiveRecord::Migration[7.0]
  def change
    add_column :jobs, :user, :string
    add_column :jobs, :start_time, :datetime
    add_column :jobs, :end_time, :datetime
  end
end
