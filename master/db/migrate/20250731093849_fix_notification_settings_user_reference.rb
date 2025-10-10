class FixNotificationSettingsUserReference < ActiveRecord::Migration[7.0]
  def up
    # Drop the existing table if it exists
    drop_table :notification_settings if table_exists?(:notification_settings)
    
    # Recreate the table with correct column types
    create_table :notification_settings, id: :integer do |t|
      t.integer :user_id, null: false
      t.boolean :notification_enabled
      t.datetime :last_notification_date
      t.datetime :last_error_date
      t.datetime :last_warning_date

      t.timestamps
    end
    
    add_index :notification_settings, :user_id
    add_foreign_key :notification_settings, :users
  end

  def down
    drop_table :notification_settings
  end
end
