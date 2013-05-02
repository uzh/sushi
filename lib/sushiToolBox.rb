#!/usr/bin/env ruby
# encoding: utf-8
module SushiToolBox
Version = '20130502-151016'

require 'pp'
require 'csv'


#require "active_record"
#  SUSHI_APP_DIR='/srv/GT/analysis/masaomi/sushi/work_party'

#  ActiveRecord::Base.establish_connection(
#              :adapter  => 'sqlite3',
#              :database => "#{SUSHI_APP_DIR}/db/development.sqlite3" 
#          )

#  require "#{SUSHI_APP_DIR}/app/models/project"
#  require "#{SUSHI_APP_DIR}/app/models/data_set"
#  require "#{SUSHI_APP_DIR}/app/models/sample"

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

include SushiToolBox

if __FILE__ == $0

  unless input_tsv=ARGV[0] and input_tsv=~/\.tsv/
    puts "Usage:\n ruby #{__FILE__} [input.tsv]"
    exit
  end

  csv = CSV.readlines(input_tsv, :col_sep=>"\t")
  data_set = []
  headers = []
  rows = []
  csv.each do |row|
    if data_set.empty?
      data_set = row
    elsif headers.empty?
      headers = row
    elsif !row.empty?
      rows << row    
    else
      save_data_set(data_set, headers, rows)
    end

    if row.empty?
      data_set = []
      headers = []
      rows = []
    end
  end

end
