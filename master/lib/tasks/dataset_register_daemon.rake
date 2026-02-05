# frozen_string_literal: true

# Rake tasks for managing the Dataset Registration Daemon
#
# Usage:
#   bundle exec rake daemon:status
#   bundle exec rake daemon:start
#   bundle exec rake daemon:stop
#   bundle exec rake daemon:restart

namespace :daemon do
  desc "Check the status of the Dataset Registration Daemon"
  task status: :environment do
    unless defined?(DatasetRegisterDaemon)
      puts "DatasetRegisterDaemon is only available in production environment"
      puts "Current environment: #{Rails.env}"
      exit 1
    end

    status = DatasetRegisterDaemon.status
    if status[:running]
      puts "Status: RUNNING"
      puts "PID: #{status[:pid]}"
    else
      puts "Status: NOT RUNNING"
      puts "Reason: #{status[:reason]}"
      puts "Stale PID: #{status[:stale_pid]}" if status[:stale_pid]
    end
  end

  desc "Start the Dataset Registration Daemon"
  task start: :environment do
    unless defined?(DatasetRegisterDaemon)
      puts "DatasetRegisterDaemon is only available in production environment"
      puts "Current environment: #{Rails.env}"
      exit 1
    end

    result = DatasetRegisterDaemon.start
    case result
    when :started
      status = DatasetRegisterDaemon.status
      puts "Daemon started successfully"
      puts "PID: #{status[:pid]}"
    when :already_running
      status = DatasetRegisterDaemon.status
      puts "Daemon is already running"
      puts "PID: #{status[:pid]}"
    when :already_locked
      puts "Another process is currently starting the daemon"
    else
      puts "Unknown result: #{result}"
    end
  end

  desc "Stop the Dataset Registration Daemon"
  task stop: :environment do
    unless defined?(DatasetRegisterDaemon)
      puts "DatasetRegisterDaemon is only available in production environment"
      puts "Current environment: #{Rails.env}"
      exit 1
    end

    result = DatasetRegisterDaemon.stop
    case result
    when :stopped
      puts "Daemon stopped successfully"
    when :not_running
      puts "Daemon was not running"
    else
      puts "Unknown result: #{result}"
    end
  end

  desc "Restart the Dataset Registration Daemon"
  task restart: :environment do
    unless defined?(DatasetRegisterDaemon)
      puts "DatasetRegisterDaemon is only available in production environment"
      puts "Current environment: #{Rails.env}"
      exit 1
    end

    puts "Stopping daemon..."
    DatasetRegisterDaemon.stop
    sleep 2

    puts "Starting daemon..."
    result = DatasetRegisterDaemon.start
    case result
    when :started
      status = DatasetRegisterDaemon.status
      puts "Daemon restarted successfully"
      puts "PID: #{status[:pid]}"
    else
      puts "Failed to start daemon: #{result}"
    end
  end
end
