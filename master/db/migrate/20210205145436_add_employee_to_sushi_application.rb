class AddEmployeeToSushiApplication < ActiveRecord::Migration[5.2]
  def change
    add_column :sushi_applications, :employee, :boolean
  end
end
