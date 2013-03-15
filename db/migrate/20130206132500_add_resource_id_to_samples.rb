class AddResourceIdToSamples < ActiveRecord::Migration
  def self.up
    add_column :samples, :resource_id, :string, :null => false, :default => "0"
  end

  def self.down
    remove_column :samples, :resource_id
  end
end
