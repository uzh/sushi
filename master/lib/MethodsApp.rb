# Methods auto-generation job.
#
# Generates an independent Slurm job that calls EzRun's generate_methods() to produce
# a Methods document after a completed analysis run.
#
# ⚠ NOT exposed in the SUSHI UI — explicitly excluded from all_sushi_applications
#   in application_controller.rb. When EzRun's generate_methods() is integrated across
#   all app classes, remove 'MethodsApp.rb' from the non_sushi_apps list there to
#   promote this to a user-submittable app.

class MethodsApp < SushiFabric::SushiApp
  def initialize(ezrun_class_name:, next_dataset_id:, gstore_result_dir:,
                 scratch_result_dir:, job_script_dir:, gstore_script_dir:,
                 sushi_server:, logger: nil, user: nil)
    super()
    @ezrun_class_name        = ezrun_class_name
    @next_dataset_id         = next_dataset_id
    @gstore_result_dir       = gstore_result_dir
    @scratch_result_dir      = scratch_result_dir
    @job_script_dir          = job_script_dir
    @gstore_script_dir       = gstore_script_dir
    @sushi_server            = sushi_server
    @logger                  = logger
    @user                    = user
    @params['process_mode']  = 'DATASET'
    @modules                 = ['AI/llm_methods_caller']
    @last_job                = true
    @queue                   = 'light'
  end

  def next_dataset
    {}
  end

  def commands
    command = ''
    command << "export LANG=en_US.UTF-8\n"
    command << "export LC_ALL=en_US.UTF-8\n"
    command << "R --vanilla --slave<<  EOT\n"
    command << "# DEV ONLY: remove .libPaths override before production\n"
    command << ".libPaths(c('/home/rdomi/R/libs', .libPaths()))\n"
    command << "if (!library(ezRun, logical.return = TRUE)){\n"
    command << "  message('retry loading ezRun')\n"
    command << "  Sys.sleep(120)\n"
    command << "  library(ezRun)\n"
    command << "}\n"
    command << "#{@ezrun_class_name}\\$new()\\$generate_methods(\n"
    command << "  gstore_script_dir = '#{@gstore_script_dir}',\n"
    command << "  output_dir        = '${SCRATCH_DIR}'\n"
    command << ")\n"
    command << "EOT\n"
    command
  end

  def job_footer
    @out.print "#### METHODS JOB DONE - COPY OUTPUTS TO GSTORE AND CLEAN UP\n"
    md_file = "methods.txt"
    @out.print copy_commands(md_file, @gstore_result_dir, 'now').join("\n"), "\n"
    @out.print <<-EOF
cd #{SushiFabric::SCRATCH_DIR}
rm -rf #{@scratch_dir} || exit 1

    EOF
  end

  def generate_script
    @job_script = File.join(@job_script_dir, "methods_dataset_#{@next_dataset_id}.sh")
    make_job_script
    @job_script
  end
end
