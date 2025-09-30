class FixNotificationSettingsUserReference < ActiveRecord::Migration[7.0]
  def up
    if table_exists?(:notification_settings)
      # Existing table: fix the column type while preserving data
      remove_foreign_key :notification_settings, :users if foreign_key_exists?(:notification_settings, :users)
      remove_index :notification_settings, :user_id if index_exists?(:notification_settings, :user_id)
      
      change_column :notification_settings, :user_id, :bigint, null: false
      
      add_index :notification_settings, :user_id
      add_foreign_key :notification_settings, :users
    else
      # New table: create with correct column types
      create_table :notification_settings do |t|
        t.bigint :user_id, null: false
        t.boolean :notification_enabled
        t.datetime :last_notification_date
        t.datetime :last_error_date
        t.datetime :last_warning_date

        t.timestamps
      end
      
      add_index :notification_settings, :user_id
      add_foreign_key :notification_settings, :users
    end
  end

  def down
    if table_exists?(:notification_settings)
      # Revert column type change if needed
      remove_foreign_key :notification_settings, :users if foreign_key_exists?(:notification_settings, :users)
      remove_index :notification_settings, :user_id if index_exists?(:notification_settings, :user_id)
      
      change_column :notification_settings, :user_id, :integer
      
      add_index :notification_settings, :user_id
    end
  end
end
