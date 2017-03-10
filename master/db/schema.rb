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
# It's strongly recommended that you check this file into your version control system.

ActiveRecord::Schema.define(version: 20170117102820) do

  create_table "data_sets", force: :cascade do |t|
    t.integer  "project_id",        limit: 4
    t.integer  "parent_id",         limit: 4
    t.string   "name",              limit: 255
    t.string   "md5",               limit: 255
    t.datetime "created_at",                                      null: false
    t.datetime "updated_at",                                      null: false
    t.string   "comment",           limit: 255
    t.text     "runnable_apps",     limit: 65535
    t.boolean  "refreshed_apps"
    t.integer  "num_samples",       limit: 4
    t.integer  "completed_samples", limit: 4
    t.integer  "user_id",           limit: 4
    t.boolean  "child",                           default: false, null: false
    t.integer  "bfabric_id",        limit: 4
  end

  create_table "jobs", force: :cascade do |t|
    t.integer  "submit_job_id",    limit: 4
    t.integer  "input_dataset_id", limit: 4
    t.integer  "next_dataset_id",  limit: 4
    t.datetime "created_at",                 null: false
    t.datetime "updated_at",                 null: false
  end

  create_table "projects", force: :cascade do |t|
    t.integer  "number",        limit: 4
    t.datetime "created_at",                     null: false
    t.datetime "updated_at",                     null: false
    t.text     "data_set_tree", limit: 16777215
  end

  create_table "samples", force: :cascade do |t|
    t.text     "key_value",   limit: 65535
    t.integer  "data_set_id", limit: 4
    t.datetime "created_at",                null: false
    t.datetime "updated_at",                null: false
  end

  create_table "sushi_applications", force: :cascade do |t|
    t.string   "class_name",        limit: 255
    t.string   "analysis_category", limit: 255
    t.text     "required_columns",  limit: 65535
    t.text     "next_dataset_keys", limit: 65535
    t.datetime "created_at",                      null: false
    t.datetime "updated_at",                      null: false
    t.text     "description",       limit: 65535
  end

  create_table "users", force: :cascade do |t|
    t.integer  "sign_in_count",       limit: 4,   default: 0
    t.datetime "current_sign_in_at"
    t.datetime "last_sign_in_at"
    t.string   "current_sign_in_ip",  limit: 255
    t.string   "last_sign_in_ip",     limit: 255
    t.integer  "selected_project",    limit: 4,   default: -1
    t.datetime "created_at",                                   null: false
    t.datetime "updated_at",                                   null: false
    t.datetime "remember_created_at"
    t.string   "login",               limit: 255, default: "", null: false
  end

  add_index "users", ["login"], name: "index_users_on_login", unique: true, using: :btree

end
