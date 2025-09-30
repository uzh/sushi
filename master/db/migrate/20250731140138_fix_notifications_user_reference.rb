class FixNotificationsUserReference < ActiveRecord::Migration[7.0]
  def up
    if table_exists?(:notifications)
      # Existing table: adjust column type safely while preserving data
      remove_foreign_key :notifications, :users if foreign_key_exists?(:notifications, :users)
      remove_index :notifications, :user_id if index_exists?(:notifications, :user_id)

      change_column :notifications, :user_id, :bigint, null: false

      add_index :notifications, :user_id
      add_foreign_key :notifications, :users
    else
      # New table: create with correct column types
      create_table :notifications do |t|
        t.bigint :user_id, null: false
        t.text :message
        t.string :notification_type
        t.boolean :read

        t.timestamps
      end

      add_index :notifications, :user_id
      add_foreign_key :notifications, :users
    end
  end

  def down
    if table_exists?(:notifications)
      remove_foreign_key :notifications, :users if foreign_key_exists?(:notifications, :users)
      remove_index :notifications, :user_id if index_exists?(:notifications, :user_id)

      change_column :notifications, :user_id, :integer

      add_index :notifications, :user_id
    end
  end
end
