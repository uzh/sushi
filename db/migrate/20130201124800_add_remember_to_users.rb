class AddRememberToUsers < ActiveRecord::Migration
  def self.up
    add_column :users, :remember_created_at, :datetime, :null => false, :default => "0000-00-00"
  end

  def self.down
    remove_column :users, :remember_created_at
  end
end
