# encoding: UTF-8
# This file is auto-generated from the current state of the database. Instead
# of editing this file, please use the migrations feature of Active Record to
# incrementally modify your database, and then regenerate this schema definition.
#
# Note that this schema.rb definition is the authoritative source for your
# database schema. If you need to create the application database on another
# system, you should be using db:schema:load, not running all the migrations
# from scratch. The latter is a flawed and unsustainable approach (the more migrations
# you'll amass, the slower it'll run and the greater likelihood for issues).
#
# It's strongly recommended to check this file into your version control system.

ActiveRecord::Schema.define(:version => 20150521082651) do

  create_table "data_sets", :force => true do |t|
    t.integer  "project_id"
    t.integer  "parent_id"
    t.string   "name"
    t.string   "md5"
    t.datetime "created_at",     :null => false
    t.datetime "updated_at",     :null => false
    t.string   "comment"
    t.text     "runnable_apps"
    t.boolean  "refreshed_apps"
  end

  create_table "jobs", :force => true do |t|
    t.integer  "submit_job_id"
    t.integer  "input_dataset_id"
    t.integer  "next_dataset_id"
    t.datetime "created_at",       :null => false
    t.datetime "updated_at",       :null => false
  end

  create_table "projects", :force => true do |t|
    t.integer  "number"
    t.datetime "created_at", :null => false
    t.datetime "updated_at", :null => false
  end

  create_table "samples", :force => true do |t|
    t.string   "key_value"
    t.integer  "data_set_id"
    t.datetime "created_at",  :null => false
    t.datetime "updated_at",  :null => false
  end

  create_table "sushi_applications", :force => true do |t|
    t.string   "class_name"
    t.string   "analysis_category"
    t.text     "required_columns"
    t.text     "next_dataset_keys"
    t.datetime "created_at",        :null => false
    t.datetime "updated_at",        :null => false
  end

  create_table "users", :force => true do |t|
    t.integer  "sign_in_count",       :default => 0
    t.datetime "current_sign_in_at"
    t.datetime "last_sign_in_at"
    t.string   "current_sign_in_ip"
    t.string   "last_sign_in_ip"
    t.integer  "selected_project",    :default => -1
    t.datetime "created_at",                          :null => false
    t.datetime "updated_at",                          :null => false
    t.datetime "remember_created_at"
    t.string   "login",               :default => "", :null => false
  end

  add_index "users", ["login"], :name => "index_users_on_login", :unique => true

end
