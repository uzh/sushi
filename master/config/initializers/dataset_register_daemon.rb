# frozen_string_literal: true

# Dataset Registration Daemon - Auto-starts on Rails boot (production only)
# Spawns lib/tasks/run_dataset_register.rb --daemon-mode as a subprocess
# Uses flock for singleton guarantee across multiple Passenger workers

if Rails.env.production?
  class DatasetRegisterDaemon
    PID_FILE = Rails.root.join('tmp', 'pids', 'dataset_register_daemon.pid')
    LOCK_FILE = Rails.root.join('tmp', 'pids', 'dataset_register_daemon.lock')

    class << self
      # Start the daemon (singleton guaranteed by flock)
      # Returns :started, :already_locked, or :already_running
      def start
        FileUtils.mkdir_p(Rails.root.join('tmp', 'pids'))

        # Try to acquire exclusive lock (singleton guarantee)
        lock_file = File.open(LOCK_FILE, File::RDWR | File::CREAT, 0644)

        unless lock_file.flock(File::LOCK_EX | File::LOCK_NB)
          # Another worker process already has the lock
          lock_file.close
          return :already_locked
        end

        begin
          # Check if daemon already running
          if daemon_already_running?
            lock_file.flock(File::LOCK_UN)
            lock_file.close
            return :already_running
          end

          # Start the daemon subprocess
          pid = spawn_daemon

          # Write PID file (safe because we hold the lock)
          File.write(PID_FILE, pid.to_s)

          Rails.logger.info "[DatasetRegisterDaemon] Started with PID #{pid}"

          # Keep lock file open to prevent other workers from starting daemon
          # Store in class variable to prevent GC from closing it
          @lock_file = lock_file

          :started
        rescue => e
          lock_file.flock(File::LOCK_UN)
          lock_file.close
          raise e
        end
      end

      # Stop the daemon
      # Returns :stopped or :not_running
      def stop
        return :not_running unless File.exist?(PID_FILE)

        pid = File.read(PID_FILE).to_i
        if process_running?(pid)
          Process.kill('TERM', pid)
          sleep 1
          Process.kill('KILL', pid) if process_running?(pid)
          Rails.logger.info "[DatasetRegisterDaemon] Stopped PID #{pid}"
        end

        FileUtils.rm_f(PID_FILE)
        :stopped
      end

      # Check daemon status
      # Returns hash with :running, :pid, :reason keys
      def status
        return { running: false, reason: 'no_pid_file' } unless File.exist?(PID_FILE)

        pid = File.read(PID_FILE).to_i
        if process_running?(pid)
          { running: true, pid: pid }
        else
          { running: false, reason: 'process_not_found', stale_pid: pid }
        end
      end

      private

      def daemon_already_running?
        return false unless File.exist?(PID_FILE)

        pid = File.read(PID_FILE).to_i
        process_running?(pid)
      end

      def process_running?(pid)
        return false if pid <= 0
        Process.kill(0, pid) # Signal 0 = check existence only
        true
      rescue Errno::ESRCH # No such process
        false
      rescue Errno::EPERM # Permission denied (but process exists)
        true
      end

      def spawn_daemon
        year = Time.now.strftime("%Y")
        script = Rails.root.join('lib', 'tasks', 'run_dataset_register.rb')
        log_file = Rails.root.join('log', "run_dataset_register_#{year}.log")
        err_file = Rails.root.join('log', "run_dataset_register_#{year}.err")

        pid = Process.spawn(
          { 'RAILS_ENV' => 'production' },
          'bundle', 'exec', 'ruby', script.to_s,
          '--run', '--daemon-mode', '--rails-root', Rails.root.to_s,
          out: [log_file.to_s, 'a'],
          err: [err_file.to_s, 'a'],
          chdir: Rails.root.to_s,
          pgroup: true # Create new process group (detach from parent)
        )

        Process.detach(pid)
        pid
      end
    end
  end

  # Start daemon after Rails fully initializes
  Rails.application.config.after_initialize do
    Thread.new do
      sleep 5 # Wait for Rails to fully initialize
      result = DatasetRegisterDaemon.start
      Rails.logger.info "[DatasetRegisterDaemon] Init result: #{result}"
    rescue => e
      Rails.logger.error "[DatasetRegisterDaemon] Failed to start: #{e.message}"
      Rails.logger.error e.backtrace.first(5).join("\n")
    end
  end
end
