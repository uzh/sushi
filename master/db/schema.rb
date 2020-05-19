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

ActiveRecord::Schema.define(version: 2017_11_24_130223) do

  create_table "data_sets", force: :cascade do |t|
    t.integer "project_id"
    t.integer "parent_id"
    t.string "name"
    t.string "md5"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.string "comment"
    t.text "runnable_apps"
    t.boolean "refreshed_apps"
    t.integer "num_samples"
    t.integer "completed_samples"
    t.integer "user_id"
    t.boolean "child", default: false, null: false
    t.integer "bfabric_id"
    t.string "sushi_app_name"
    t.string "run_name_order_id"
    t.integer "workunit_id"
  end

  create_table "jobs", force: :cascade do |t|
    t.integer "submit_job_id"
    t.integer "input_dataset_id"
    t.integer "next_dataset_id"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.string "script_path"
  end

  create_table "projects", force: :cascade do |t|
    t.integer "number"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.text "data_set_tree", limit: 16777215
  end

  create_table "samples", force: :cascade do |t|
    t.text "key_value"
    t.integer "data_set_id"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
  end

  create_table "sushi_applications", force: :cascade do |t|
    t.string "class_name"
    t.string "analysis_category"
    t.text "required_columns"
    t.text "next_dataset_keys"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.text "description"
  end

  create_table "users", force: :cascade do |t|
    t.integer "sign_in_count", default: 0
    t.datetime "current_sign_in_at"
    t.datetime "last_sign_in_at"
    t.string "current_sign_in_ip"
    t.string "last_sign_in_ip"
    t.integer "selected_project", default: -1
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.datetime "remember_created_at"
    t.string "login", default: "", null: false
    t.index ["login"], name: "index_users_on_login", unique: true
  end

end
