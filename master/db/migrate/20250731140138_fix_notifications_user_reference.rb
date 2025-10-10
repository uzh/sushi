class FixNotificationsUserReference < ActiveRecord::Migration[7.0]
  def up
    # Drop the existing table if it exists
    drop_table :notifications if table_exists?(:notifications)
    
    # Recreate the table with correct column types
    create_table :notifications, id: :integer do |t|
      t.integer :user_id, null: false
      t.text :message
      t.string :notification_type
      t.boolean :read

      t.timestamps
    end
    
    add_index :notifications, :user_id
    add_foreign_key :notifications, :users
  end

  def down
    drop_table :notifications
  end
end
