class AddSubmitCommandAndStatusToJobs < ActiveRecord::Migration[7.0]
  def change
    add_column :jobs, :submit_command, :text
    add_column :jobs, :status, :string
  end
end
