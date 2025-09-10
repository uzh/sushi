#!/usr/bin/env ruby
# encoding: utf-8
# Version = '20250730-053958'

require 'pathname'
require 'active_support/core_ext/time'

RAILS_ROOT='/srv/sushi/production/master'
YEAR = Time.now.strftime("%Y")

# Load Rails environment for notification service
def load_rails_environment
  root = (defined?($RUNNER_RAILS_ROOT) && $RUNNER_RAILS_ROOT && !$RUNNER_RAILS_ROOT.to_s.empty?) ? $RUNNER_RAILS_ROOT : RAILS_ROOT
  require 'bundler/setup'
  require File.join(root, 'config', 'environment')
rescue LoadError => e
  puts "Warning: Could not load Rails environment: #{e.message}"
  puts "Notifications will not be sent."
  return false
end

# Send notification for errors or warnings
def send_notification(error_message = nil, warning_message = nil)
  puts "Sending notification - Error: #{error_message}, Warning: #{warning_message}"
  return unless load_rails_environment
  
  begin
    notification_service = NotificationService.new
    notification_service.notify_dataset_registration_issues(error_message, warning_message)
    puts "Notification sent successfully"
  rescue => e
    puts "Error sending notification: #{e.message}"
  end
end

help =-> () do
  puts <<-eos
  usage:
   #{File.basename(__FILE__)} [task name] (--run, --year [YEAR], --rails-root [RAILS_ROOT], --daemon-mode, --help)

  RAILS_ROOT (default): #{RAILS_ROOT}

  default (for test run):
   #{File.basename(__FILE__)} --year #{YEAR}

  recommend (for production run):
   #{File.basename(__FILE__)} --run --daemon-mode > #{File.basename(__FILE__).gsub(/.rb/, "_#{Time.now.strftime("%Y%m%d")}.log")}

  --daemon-mode
   * start the registraiton process at 0:00 AM every night

  e.g.:
   #{File.basename(__FILE__)} --year 2024
   #{File.basename(__FILE__)} --run
   #{File.basename(__FILE__)} --year 2024 --run
   #{File.basename(__FILE__)} ds:register_datasets[2024,]
   #{File.basename(__FILE__)} ds:register_datasets[2024,run]
  eos
  exit
end

if ARGV.index("-h") or ARGV.index("--help")
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
unless File.exist?(rails_root)
  warn "# WARN: #{rails_root} does not exist"
  help.()
end

# share rails_root with notification loader
$RUNNER_RAILS_ROOT = rails_root.to_s

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

# Run command and capture stdout/stderr separately; notify if stderr seen or non-zero exit
def run_and_monitor(command)
  r_out, w_out = IO.pipe
  r_err, w_err = IO.pipe
  pid = Process.spawn(command, out: w_out, err: w_err)
  w_out.close
  w_err.close

  stderr_buf = +""
  t_out = Thread.new { r_out.each_line { |l| puts l } }
  t_err = Thread.new { r_err.each_line { |l| STDERR.puts l; stderr_buf << l } }

  _, status = Process.wait2(pid)
  t_out.join
  t_err.join
  [status.success?, stderr_buf, status.exitstatus]
end
Dir.chdir(rails_root) do
  if daemon_mode
    loop do
      begin
        now = Time.now
        # Calculate next midnight using a safer method
        next_midnight = (now + 86400).beginning_of_day
        wait_time = next_midnight - now
        sleep(wait_time)
        
        # Execute the command and capture output
        puts "Executing dataset registration task at #{Time.now}"
        success, stderr_buf, exitstatus = run_and_monitor(command)
        if !success || !stderr_buf.strip.empty?
          msg = "Dataset registration task issue: exit=#{exitstatus}, stderr:\n" + stderr_buf.lines.first(200).join
          puts msg
          send_notification(msg)
        else
          puts "Dataset registration task completed successfully"
        end
      rescue ArgumentError => e
        # Fallback when Time.new argument error occurs
        puts "Warning: Time calculation error occurred: #{e.message}"
        puts "Using fallback method to calculate next midnight..."
        now = Time.now
        # Move 24 hours from current time and then set to 0:00
        next_midnight = Time.at(now.to_i + 86400).beginning_of_day
        wait_time = next_midnight - now
        sleep(wait_time)
        
        puts "Executing dataset registration task at #{Time.now}"
        success, stderr_buf, exitstatus = run_and_monitor(command)
        if !success || !stderr_buf.strip.empty?
          msg = "Dataset registration task issue: exit=#{exitstatus}, stderr:\n" + stderr_buf.lines.first(200).join
          puts msg
          send_notification(msg)
        else
          puts "Dataset registration task completed successfully"
        end
      rescue => e
        error_message = "Error in daemon mode: #{e.message}"
        puts error_message
        send_notification(error_message)
        puts "Waiting 60 seconds before retrying..."
        sleep(60)
      end
    end
  else
    # Execute the command and capture output
    puts "Executing dataset registration task at #{Time.now}"
    success, stderr_buf, exitstatus = run_and_monitor(command)
    if !success || !stderr_buf.strip.empty?
      msg = "Dataset registration task issue: exit=#{exitstatus}, stderr:\n" + stderr_buf.lines.first(200).join
      puts msg
      send_notification(msg)
    else
      puts "Dataset registration task completed successfully"
    end
  end
end


