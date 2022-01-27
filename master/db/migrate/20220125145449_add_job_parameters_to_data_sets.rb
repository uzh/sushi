class AddJobParametersToDataSets < ActiveRecord::Migration[5.2]
  def change
    add_column :data_sets, :job_parameters, :text
  end
end
