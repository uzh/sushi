# Methods auto-generation job.
#
# Generates an independent Slurm job that calls EzRun's write_methods() to produce
# a Methods document after a completed analysis run.
#
# ⚠ NOT exposed in the SUSHI UI — explicitly excluded from all_sushi_applications
#   in application_controller.rb. When EzRun's write_methods() is integrated across
#   all app classes, remove 'MethodsApp.rb' from the non_sushi_apps list there to
#   promote this to a user-submittable app.

require 'shellwords'

class MethodsApp < SushiFabric::SushiApp
  # Set by hand -- update if the model backing AI/llm_methods_caller changes.
  # Source: LLM_CALLER_MODEL in the AI/llm_methods_caller module, which llm_write_methods
  # uses as its --model default when SUSHI doesn't override it (it never does).
  LLM_MODEL_NAME = 'DeepSeek-V4-Flash-DSpark'

  def initialize(ezrun_class_name:, analysis_name:, next_dataset_id:, gstore_result_dir:,
                 scratch_result_dir:, job_script_dir:, gstore_script_dir:,
                 sushi_server:, logger: nil, user: nil, parent_methods_path: nil, citation: [])
    super()
    @ezrun_class_name        = ezrun_class_name
    @analysis_name           = analysis_name
    @citation                = citation
    @next_dataset_id         = next_dataset_id
    @gstore_result_dir       = gstore_result_dir
    @scratch_result_dir      = scratch_result_dir
    @job_script_dir          = job_script_dir
    @gstore_script_dir       = gstore_script_dir
    @sushi_server            = sushi_server
    @logger                  = logger
    @user                    = user
    @parent_methods_path     = parent_methods_path
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

    has_citations = @citation && !@citation.empty?
    citations_file = '${SCRATCH_DIR}/citations_candidates.txt'
    if has_citations
      # A plain local file, written before R starts, keeps the actual bibliographic
      # text (accents, ampersands, parens) out of the R heredoc entirely -- same
      # Shellwords-escaped-printf technique already used for references_text below.
      command << "printf '%s\\n' #{Shellwords.escape(@citation.join("\n"))} > #{citations_file}\n"
    end

    command << "R --vanilla --slave<<  EOT\n"
    command << "if (!library(ezRun, logical.return = TRUE)){\n"
    command << "  message('retry loading ezRun')\n"
    command << "  Sys.sleep(120)\n"
    command << "  library(ezRun)\n"
    command << "}\n"
    command << "#{@ezrun_class_name}\\$new()\\$write_methods(\n"
    args = ["  gstore_script_dir = '#{@gstore_script_dir}'", "  output_dir        = '${SCRATCH_DIR}'"]
    args << "  citations         = readLines('#{citations_file}')" if has_citations
    command << args.join(",\n") << "\n"
    command << ")\n"
    command << "EOT\n"
    command
  end

  def job_footer
    @out.print "#### METHODS JOB DONE - COPY OUTPUTS TO GSTORE AND CLEAN UP\n"
    md_file = "methods.md"

    # write_methods() (and any app-specific override) only ever writes its own prose -- it
    # has no idea about ancestry, its app name, or the references/declaration boilerplate.
    # Wrap and chain the section here, once, so no per-app override can forget any of it.
    if @citation && !@citation.empty?
      # The model was asked (llm_write_methods' REFERENCES_INSTRUCTION) to end its
      # response with a "## References" line, then a verbatim copy of whichever
      # candidates it found actual evidence for. Never trust that copy as the source
      # of truth: for each known citation, check whether its DOI/URL -- the span an
      # LLM is least likely to mangle -- appears anywhere in the raw response. A hit
      # emits OUR OWN stored string, never the model's; a miss drops it. Marker
      # absent entirely (old ezRun, or the model didn't comply) falls back to the
      # full static list -- never worse than before this feature existed.
      anchors        = @citation.map { |c| c[%r{https?://\S+}] || c }
      citations_bash = @citation.map { |c| Shellwords.escape(c) }.join(' ')
      anchors_bash   = anchors.map { |a| Shellwords.escape(a) }.join(' ')
      full_list_bash = Shellwords.escape(@citation.join("\n"))
      @out.print <<-EOF
body="$(cat #{md_file})"
if grep -qF -- '## References' <<< "$body"; then
  description_body="$(awk '/^## References/{exit} {print}' <<< "$body")"
  citations=(#{citations_bash})
  anchors=(#{anchors_bash})
  kept=()
  for i in "${!citations[@]}"; do
    if grep -qF -- "${anchors[$i]}" <<< "$body"; then
      kept+=("${citations[$i]}")
    fi
  done
  if [ ${#kept[@]} -gt 0 ]; then
    references_text="$(printf '%s\\n' "${kept[@]}")"
  else
    references_text="pending"
  fi
else
  description_body="$body"
  references_text=#{full_list_bash}
fi
      EOF
    else
      @out.print <<-EOF
body="$(cat #{md_file})"
description_body="$body"
references_text="pending"
      EOF
    end
    @out.print <<-EOF
{
  printf '## #{@analysis_name} | %s\\n\\n' "$(date '+%Y-%m-%d %H:%M')"
  printf '### Description\\n\\n'
  printf '%s\\n\\n' "$description_body"
  printf '### References\\n\\n'
  printf '%s\\n\\n' "$references_text"
  printf '### Declaration\\n\\n'
  printf 'This description has been generated by #{LLM_MODEL_NAME} based on the input data, the parameters and the result of the analysis.\\n'
} > #{md_file}.section
mv #{md_file}.section #{md_file}
    EOF

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
    @job_script = File.join(@job_script_dir, "methods_dataset_#{@next_dataset_id}.sh")
    make_job_script
    @job_script
  end
end
