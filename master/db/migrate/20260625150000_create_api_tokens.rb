class CreateApiTokens < ActiveRecord::Migration[7.0]
  def change
    create_table :api_tokens do |t|
      t.string   :token_hash, null: false
      t.string   :name
      t.text     :scope         # serialized array of project numbers
      t.datetime :expires_at
      t.datetime :revoked_at
      t.timestamps
    end
    add_index :api_tokens, :token_hash, unique: true
  end
end
