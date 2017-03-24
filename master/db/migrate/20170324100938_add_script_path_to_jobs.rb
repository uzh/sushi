class AddScriptPathToJobs < ActiveRecord::Migration
  def change
    add_column :jobs, :script_path, :string
  end
end
