# Methods auto-generation job.
#
# Generates an independent Slurm job that calls EzRun's write_methods() to produce
# a Methods document after a completed analysis run.
#
# ⚠ NOT exposed in the SUSHI UI — explicitly excluded from all_sushi_applications
#   in application_controller.rb. When EzRun's write_methods() is integrated across
#   all app classes, remove 'MethodsApp.rb' from the non_sushi_apps list there to
#   promote this to a user-submittable app.

class MethodsApp < SushiFabric::SushiApp
  # Per-app documentation lives next to the app classes, one file per SUSHI app.
  APP_DOCS_DIR = File.join(__dir__, 'docs')

  def initialize(ezrun_class_name:, analysis_name:, next_dataset_id:, gstore_result_dir:,
                 scratch_result_dir:, job_script_dir:, gstore_script_dir:,
                 sushi_server:, logger: nil, user: nil, parent_methods_path: nil,
                 sample_count: nil, example_script: nil, sushi_app_name: nil)
    super()
    @ezrun_class_name        = ezrun_class_name
    @sushi_app_name          = sushi_app_name
    @analysis_name           = analysis_name
    @next_dataset_id         = next_dataset_id
    @gstore_result_dir       = gstore_result_dir
    @scratch_result_dir      = scratch_result_dir
    @job_script_dir          = job_script_dir
    @gstore_script_dir       = gstore_script_dir
    @sushi_server            = sushi_server
    @logger                  = logger
    @user                    = user
    @parent_methods_path     = parent_methods_path
    @sample_count            = sample_count
    @example_script          = example_script
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
    command << "if (!library(ezRun, logical.return = TRUE)){\n"
    command << "  message('retry loading ezRun')\n"
    command << "  Sys.sleep(120)\n"
    command << "  library(ezRun)\n"
    command << "}\n"
    command << "#{@ezrun_class_name}\\$new()\\$write_methods(\n"
    args = ["  gstore_script_dir = '#{@gstore_script_dir}'",
            "  output_dir        = '${SCRATCH_DIR}'",
            "  analysis_name     = '#{@analysis_name}'"]
    args << "  example_script    = '#{@example_script}'" if @example_script
    args << "  sample_count      = #{@sample_count}" if @sample_count
    args << "  app_doc           = '#{@app_doc_path}'" if @app_doc_path
    command << args.join(",\n") << "\n"
    command << ")\n"
    command << "EOT\n"
    command
  end

  def job_footer
    @out.print "#### METHODS JOB DONE - COPY OUTPUTS TO GSTORE AND CLEAN UP\n"
    md_file = "methods.md"

    # write_methods() already wrote the complete section (heading, Description,
    # References, Declaration) -- the only thing left here is chaining onto the
    # parent's own methods.md, if one exists.
    if @parent_methods_path
      @out.print <<-EOF
if [ -f "#{@parent_methods_path}" ]; then
  { cat "#{@parent_methods_path}"; printf '\\n'; cat #{md_file}; } > #{md_file}.chained
  mv #{md_file}.chained #{md_file}
fi
      EOF
    end

    @out.print copy_commands(md_file, @gstore_result_dir, 'now').join("\n"), "\n"
    @out.print <<-EOF
cd #{SushiFabric::SCRATCH_DIR}
rm -rf #{@scratch_dir} || exit 1

    EOF
  end

  def generate_script
    stage_app_doc
    @job_script = File.join(@job_script_dir, "methods_dataset_#{@next_dataset_id}.sh")
    make_job_script
    @job_script
  end

  private

  # Copies the app's documentation next to the job scripts. It then travels to gstore
  # with them in the regular copy wave, so the methods job reads it from the run's own
  # scripts/ folder, and every run keeps the doc version its Methods text was written
  # from. An app without a doc gets no app_doc argument: the job runs as before.
  def stage_app_doc
    return unless @sushi_app_name
    doc = File.join(APP_DOCS_DIR, "#{@sushi_app_name}.md")
    return unless File.exist?(doc)
    FileUtils.cp(doc, @job_script_dir)
    @app_doc_path = File.join(@gstore_script_dir, File.basename(doc))
  rescue => e
    @logger.error("app doc staging failed: #{e.message}") if @logger
    @app_doc_path = nil
  end
end
