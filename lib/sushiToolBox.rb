#!/usr/bin/env ruby
# encoding: utf-8
Version = '20130530-141919'

require "active_record"

SUSHI_APP_DIR='/srv/GT/analysis/masaomi/sushi/work_party'
SUSHI_DB_TYPE='sqlite3'

module SushiToolBox
  begin
    ::Project
  rescue
    ActiveRecord::Base.establish_connection(
                :adapter  => SUSHI_DB_TYPE,
                :database => "#{SUSHI_APP_DIR}/db/development.sqlite3" 
            )
    require "#{SUSHI_APP_DIR}/app/models/project"
    require "#{SUSHI_APP_DIR}/app/models/data_set"
    require "#{SUSHI_APP_DIR}/app/models/sample"
  end

  def save_data_set(data_set_arr, headers, rows)
    data_set_hash = Hash[*data_set_arr]
    if project = Project.find_by_number(data_set_hash['ProjectNumber'].to_i)
      data_set = DataSet.new
      data_set.name = data_set_hash['DataSetName']
      data_set.project = project
      if parent_id = data_set_hash['ParentID'] and parent_data_set = DataSet.find_by_id(parent_id.to_i)
        data_set.data_set = parent_data_set
      end

      sample_hash = {}
      rows.each do |row|
        headers.each_with_index do |header, i|
         sample_hash[header]=row[i]
        end
        sample = Sample.new
        sample.key_value = sample_hash.to_s
        sample.save unless sample.saved?
        data_set.samples << sample
      end

      data_set.md5 = data_set.md5hexdigest
      unless data_set.saved?
        project.data_sets << data_set
        parent_data_set.data_sets << data_set if parent_data_set
        data_set.save
      end
    end
  end
end

