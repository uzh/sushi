class RemoveFinalLogPathFromJobs < ActiveRecord::Migration[7.0]
  def change
    remove_column :jobs, :final_log_path, :string
  end
end
