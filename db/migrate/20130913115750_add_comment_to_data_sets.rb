class AddCommentToDataSets < ActiveRecord::Migration
  def change
    add_column :data_sets, :comment, :string
  end
end
