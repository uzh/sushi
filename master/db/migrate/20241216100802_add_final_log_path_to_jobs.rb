class AddFinalLogPathToJobs < ActiveRecord::Migration[7.0]
  def change
    add_column :jobs, :final_log_path, :string
  end
end
