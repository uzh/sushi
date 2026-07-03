class AddPrincipalAndLoginToApiTokens < ActiveRecord::Migration[7.0]
  # Per-user (LDAP-bound) token principal. A token is either `static` (its scope
  # is the stored project-number array — existing behavior) or `user` (authorized
  # live against the bound LDAP login's FGCZ project membership). Additive only:
  # existing rows default to `static` with a NULL login, so their behavior is
  # unchanged (design v0.7 P-INV-8).
  def change
    add_column :api_tokens, :principal, :string, null: false, default: "static"
    add_column :api_tokens, :login, :string
  end
end
