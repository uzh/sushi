# encoding: utf-8
# This file should contain all the record creation needed to seed the database with its default values.
# The data can then be loaded with the rake db:seed (or created alongside the db with db:setup).
#
# Examples:
#
#   cities = City.create([{ name: 'Chicago' }, { name: 'Copenhagen' }])
#   Mayor.create(name: 'Emanuel', city: cities.first)

Sample.create(name: 'sample_name1', path:'/srv/GT/analysis/masaomi/sushi/work/public/20120905-I8_970_20081205_R1_sample.fastq.gz')
Sample.create(name: 'sample_name2', path:'/srv/GT/analysis/masaomi/sushi/work/public/20120905-I8_970_20081205_R2_sample.fastq.gz')
DataList.create(data_set_id: 1, sample_id: 1)
DataList.create(data_set_id: 2, sample_id: 2)
DataList.create(data_set_id: 1, sample_id: 2)
DataSet.create(name: 'data_set1', note: 'data set 1')
DataSet.create(name: 'data_set2', note: 'data set 2')
