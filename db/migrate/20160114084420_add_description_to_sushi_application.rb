class AddDescriptionToSushiApplication < ActiveRecord::Migration
  def change
    add_column :sushi_applications, :description, :text
  end
end
