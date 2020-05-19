class AddCommentToDataSets < ActiveRecord::Migration[5.2]
  def change
    add_column :data_sets, :comment, :string
  end
end
