#!/usr/bin/env ruby
# encoding: utf-8
# Version = '20240718-110347'

require 'pathname'

RAILS_ROOT='/srv/sushi/production/master'
SLEEP_TIME = 60*60 # [s] = 1h, valid only for daemon-mode

help =-> () do
  puts <<-eos
  usage:
   #{File.basename(__FILE__)} [task name] (--run, --year [YEAR], --rails-root [RAILS_ROOT], --daemon-mode)

  RAILS_ROOT (default): #{RAILS_ROOT}

  e.g.:
   #{File.basename(__FILE__)} --year 2024
   #{File.basename(__FILE__)} --run
   #{File.basename(__FILE__)} --year 2024 --run
   #{File.basename(__FILE__)} ds:register_datasets[2024,]
   #{File.basename(__FILE__)} ds:register_datasets[2024,run]
  eos
  exit
end

if i=ARGV.index("-h")
  help.()
end

run = ARGV.index("--run")
year = if i=ARGV.index("--year")
         ARGV[i+1]
       else
         Time.now.strftime("%Y")
       end

rails_root = if i=ARGV.index("--rails-root")
               ARGV[i+1]
             else
               Pathname.new(RAILS_ROOT)
             end

daemon_mode = ARGV.index("--daemon-mode")

ENV['RAILS_ENV'] = 'production'
ENV['DISABLE_DATABASE_ENVIRONMENT_CHECK'] = '1'

task_name_with_args = if ARGV[0] =~ /^ds\:register_datasets/
                        ARGV[0]
                      elsif run
                        "ds:register_datasets[#{year},run]"
                      else
                        "ds:register_datasets[#{year},]"
                      end

command = "bundle exec rake #{task_name_with_args}"
Dir.chdir(rails_root) do
  if daemon_mode
    loop do
      system(command)
      sleep SLEEP_TIME
    end
  else
    system(command)
  end
end


