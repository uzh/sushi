class AddDescriptionToSushiApplication < ActiveRecord::Migration[5.2]
  def change
    add_column :sushi_applications, :description, :text
  end
end
