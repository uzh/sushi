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

ActiveRecord::Schema[7.0].define(version: 2025_09_26_000000) do
  create_table "data_sets", id: :integer, charset: "utf8mb3", collation: "utf8mb3_unicode_ci", force: :cascade do |t|
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
    t.index ["parent_id"], name: "index_data_sets_on_parent_id"
    t.index ["project_id"], name: "index_data_sets_on_project_id"
  end

  create_table "jobs", id: :integer, charset: "utf8mb3", collation: "utf8mb3_unicode_ci", force: :cascade do |t|
    t.integer "submit_job_id"
    t.integer "input_dataset_id"
    t.integer "next_dataset_id"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.string "script_path"
    t.string "stdout_path"
    t.string "stderr_path"
    t.text "submit_command"
    t.string "status"
    t.string "user"
    t.datetime "start_time"
    t.datetime "end_time"
    t.index ["id", "next_dataset_id"], name: "index_jobs_on_id_and_next_dataset_id"
    t.index ["input_dataset_id"], name: "index_jobs_on_input_dataset_id"
    t.index ["next_dataset_id", "id"], name: "index_jobs_on_next_dataset_id_and_id"
    t.index ["next_dataset_id"], name: "index_jobs_on_next_dataset_id"
    t.index ["status"], name: "index_jobs_on_status"
  end

  create_table "notification_settings", charset: "utf8mb3", collation: "utf8mb3_general_ci", force: :cascade do |t|
    t.integer "user_id", null: false
    t.boolean "notification_enabled"
    t.datetime "last_notification_date"
    t.datetime "last_error_date"
    t.datetime "last_warning_date"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.index ["user_id"], name: "index_notification_settings_on_user_id"
  end

  create_table "notifications", charset: "utf8mb3", collation: "utf8mb3_general_ci", force: :cascade do |t|
    t.integer "user_id", null: false
    t.text "message"
    t.string "notification_type"
    t.boolean "read"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.index ["user_id"], name: "index_notifications_on_user_id"
  end

  create_table "projects", id: :integer, charset: "utf8mb3", collation: "utf8mb3_unicode_ci", force: :cascade do |t|
    t.integer "number"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.text "data_set_tree", size: :medium
    t.index ["number"], name: "index_projects_on_number"
  end

  create_table "samples", id: :integer, charset: "utf8mb3", collation: "utf8mb3_unicode_ci", force: :cascade do |t|
    t.text "key_value"
    t.integer "data_set_id"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.index ["data_set_id"], name: "index_samples_on_data_set_id"
  end

  create_table "sushi_applications", id: :integer, charset: "utf8mb3", collation: "utf8mb3_unicode_ci", force: :cascade do |t|
    t.string "class_name"
    t.string "analysis_category"
    t.text "required_columns"
    t.text "next_dataset_keys"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.text "description"
    t.boolean "employee"
  end

  create_table "users", id: :integer, charset: "utf8mb3", collation: "utf8mb3_unicode_ci", force: :cascade do |t|
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

  add_foreign_key "notification_settings", "users"
  add_foreign_key "notifications", "users"
end
