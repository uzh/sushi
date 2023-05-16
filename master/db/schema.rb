# This file is auto-generated from the current state of the database. Instead
# of editing this file, please use the migrations feature of Active Record to
# incrementally modify your database, and then regenerate this schema definition.
#
# This file is the source Rails uses to define your schema when running `bin/rails
# db:schema:load`. When creating a new database, `bin/rails db:schema:load` tends to
# be faster and is potentially less error prone than running all of your
# migrations from scratch. Old migrations may fail to apply correctly if those
# migrations use external dependencies or application code.
#
# It's strongly recommended that you check this file into your version control system.

ActiveRecord::Schema.define(version: 2023_05_16_090858) do

  create_table "data_sets", id: :integer, charset: "utf8", collation: "utf8_unicode_ci", force: :cascade do |t|
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
    t.text "order_ids"
    t.text "job_parameters"
    t.integer "order_id"
    t.index ["order_id"], name: "index_data_sets_on_order_id"
  end

  create_table "jobs", id: :integer, charset: "utf8", collation: "utf8_unicode_ci", force: :cascade do |t|
    t.integer "submit_job_id"
    t.integer "input_dataset_id"
    t.integer "next_dataset_id"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.string "script_path"
  end

  create_table "projects", id: :integer, charset: "utf8", collation: "utf8_unicode_ci", force: :cascade do |t|
    t.integer "number"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.text "data_set_tree", size: :medium
  end

  create_table "samples", id: :integer, charset: "utf8", collation: "utf8_unicode_ci", force: :cascade do |t|
    t.text "key_value"
    t.integer "data_set_id"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
  end

  create_table "sushi_applications", id: :integer, charset: "utf8", collation: "utf8_unicode_ci", force: :cascade do |t|
    t.string "class_name"
    t.string "analysis_category"
    t.text "required_columns"
    t.text "next_dataset_keys"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.text "description"
    t.boolean "employee"
  end

  create_table "users", id: :integer, charset: "utf8", collation: "utf8_unicode_ci", force: :cascade do |t|
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
