class AddRegistrationKeyToDataSets < ActiveRecord::Migration[7.0]
  def change
    add_column :data_sets, :registration_key, :string
    add_index  :data_sets, :registration_key, unique: true
  end
end
