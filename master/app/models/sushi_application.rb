class SushiApplication < ActiveRecord::Base
#  attr_accessible :analysis_category, :class, :next_dataset_keys, :required_columns, :description
  serialize :required_columns
  serialize :next_dataset_keys

  def required_columns_satisfied_by?(headers)
    normalized_headers = headers.map { |col| col.to_s.gsub(/\[.+\]/, '').strip }

    if required_columns.all? { |entry| entry.is_a?(Array) }
      # XOR-mode
      satisfied_options = required_columns.count do |option|
        Array(option).all? { |req| normalized_headers.include?(req.gsub(/\[.+\]/, '').strip) }
      end
      satisfied_options == 1
    else
      # AND-mode (previous)
      missing = required_columns - normalized_headers
      missing.empty?
    end
  end
end
