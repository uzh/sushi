class SushiApplication < ActiveRecord::Base
#  attr_accessible :analysis_category, :class, :next_dataset_keys, :required_columns, :description
  serialize :required_columns
  serialize :next_dataset_keys
end
