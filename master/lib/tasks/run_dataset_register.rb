#!/usr/bin/env ruby
# encoding: utf-8
# Version = '20240712-101921'

#!/usr/bin/env ruby

require 'pathname'

help =-> () do
  puts <<-eos
  usage:
   #{File.basename(__FILE__)} [task name]

  e.g.:
   #{File.basename(__FILE__)} ds:register_datasets[2024,run]
   #{File.basename(__FILE__)} ds:register_datasets[2024,]
  eos
  exit
end

if i=ARGV.index("-h")
  help.()
end


# Railsアプリケーションのディレクトリを指定（絶対パスに置き換えてください）
rails_root = Pathname.new('/srv/sushi/masa_test_sushi_20240711/master')

# 環境変数の設定
ENV['RAILS_ENV'] = 'production'
ENV['DISABLE_DATABASE_ENVIRONMENT_CHECK'] = '1'

# 実行するRakeタスクとその引数を取得
task_name_with_args = ARGV[0] || 'ds:register_datasets[2024,]'

# bundle exec rakeコマンドを実行
command = "bundle exec rake #{task_name_with_args}"
Dir.chdir(rails_root) do
  system(command)
end


# #!/usr/bin/env ruby
# 
# require 'pathname'
# require 'rake'
# 
# # Railsアプリケーションのディレクトリを指定（絶対パスに置き換えてください）
# rails_root = Pathname.new('/srv/sushi/masa_test_sushi_20240704/master')
# 
# # 環境変数の設定
# ENV['RAILS_ENV'] = 'production'
# ENV['DISABLE_DATABASE_ENVIRONMENT_CHECK'] = '1'
# 
# # Rails環境をロード
# require rails_root.join('config', 'environment')
# 
# ## Rakeアプリケーションを初期化
# #rake_app = Rake.application
# #rake_app.init
# #rake_app.load_rakefile
# 
# # 実行したいRakeタスクを指定
# # 引数からタスク名とその引数を取得
# # ds:register_datasets[2024,run]
# task_name_with_args = ARGV[0] || 'ds:register_datasets[2024,]'
# 
# 
# #task_name, args = task_name_with_args.split('[')
# #args = args ? args.chomp(']').split(',') : []
# 
# # bundle exec rakeコマンドを実行
# command = "bundle exec rake #{task_name_with_args}"
# Dir.chdir(rails_root) do
#   system(command)
# end
# 
# # # タスクを実行
# # rake_app[task_name].invoke(*args)
# # 
# # # タスクを実行
# # rake_app[task_name].invoke
# 
# 
# 
# 
