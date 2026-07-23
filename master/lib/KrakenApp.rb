#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class KrakenApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Kraken'
    @analysis_category = 'Metagenomics'
    @description =<<-EOS
Kraken taxonomic sequence classification system
<a href='https://ccb.jhu.edu/software/kraken2/index.shtml'>https://ccb.jhu.edu/software/kraken2/index.shtml</a>
<br/><br/>
By default each sample is classified in its own job (SAMPLE mode). Enable
<b>exclusive</b> to instead process the whole dataset in a single job on one node
reserved with SLURM <code>--exclusive</code>, classifying all samples sequentially
through one Kraken2 <code>--use-daemon</code> classifier so the (large) database is
loaded into memory only once for the run instead of once per sample. The output
table is the same either way (one row per sample).
EOS

    @required_columns = ['Name','Read1']
    @required_params = ['paired', 'krakenDBOpt']
    # optional params
    # 'exclusive' is the opt-in switch for the fast, database-loaded-once path.
    # OFF (default): classic SAMPLE mode — one job per sample, no node reservation.
    # ON: the whole dataset runs as a single job on a node reserved with SLURM
    # --exclusive, and the ezRun worker classifies all samples sequentially through
    # one shared k2 classify daemon (started with --use-daemon on the first sample,
    # stopped with `k2 clean --stop-daemon` after the last) so the DB loads once.
    # --exclusive is REQUIRED for that path: the daemon addresses its control channel
    # through fixed files in /tmp, which is shared between jobs on a node, so a
    # co-located second Kraken job would collide with it. process_mode is derived
    # from this flag in preprocess, and the worker only starts the daemon when it is
    # on — so the two can never be set inconsistently.
    @params['process_mode'] = 'SAMPLE'
    @params['exclusive'] = false
    @params['exclusive', 'description'] = 'load the Kraken2 database only once for the whole dataset: reserve a full node (SLURM --exclusive) and classify all samples sequentially through one shared k2 --use-daemon classifier. Off = one job per sample (classic behaviour). Recommended for large databases and/or many samples.'
    @params['exclusive', 'context'] = 'slurm'
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '100'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '500'
    @params['scratch', "context"] = "slurm"
    @params['paired'] = false
    @params['paired', "context"] = "Kraken"
    
    @params['krakenDBOpt'] = kraken_db_choices
    @params['krakenDBOpt', 'description'] = 'kraken database options with size on disk (please request at least a slightly larger amount of RAM in parameters above). Auto-detected from /srv/GT/databases/kraken2/ (see also https://benlangmead.github.io/aws-indexes/k2). With multiDB checked, Ctrl/Cmd-click to select multiple — k2 v2.17+ will classify against all in one pass (--db DB1,DB2,...). WARNING: multi-DB mode is NCBI ONLY — every selected DB must share the NCBI taxonomy tree (custom DBs like GTDB, SILVA, UHGG cannot be mixed); --report-minimizer-data is auto-disabled in multi-DB mode (kraken2 refuses it).'
    @params['krakenDBOpt', "context"] = "Kraken"
    @params['krakenDBOpt', 'multi_selection'] = true
    @params['multiDB'] = false
    @params['multiDB', 'description'] = 'enable multi-database classification (k2 v2.17, --db DB1,DB2,...). MULTI-DB IS NCBI ONLY — all selected DBs must share the NCBI taxonomy. If unchecked, only the first selected DB is used.'
    @params['multiDB', "context"] = "Kraken"
    @params['krakenConfidenceOpt'] = '0.0'
    @params['krakenConfidenceOpt', 'description'] = 'Confidence score threshold, between 0 and 1'
    @params['krakenConfidenceOpt', "context"] = "Kraken"
    @params['krakenPhredOpt'] = '0'
    @params['krakenPhredOpt', 'description'] = 'Phred score threshold'
    @params['krakenPhredOpt', "context"] = "Kraken"
    @params['minimum_hit_groups'] = ['2', '3', '4']
    @params['minimum_hit_groups', 'description'] = 'minimum number of hit groups (overlapping k-mers sharing the same minimizer) needed to make a call (kraken2 default: 2)'
    @params['minimum_hit_groups', "context"] = "Kraken"
    @params['report_minimizer_data'] = ['no', 'yes']
    @params['report_minimizer_data', 'description'] = 'add per-clade distinct-minimizer columns to the report (enables advanced filtering in exploreMetaTax)'
    @params['report_minimizer_data', "context"] = "Kraken"
    @params['save_unclassified'] = ['no', 'yes']
    @params['save_unclassified', 'description'] = 'save reads that did not classify to fastq.gz. Paired-end → &lt;sample&gt;_unclassified_{1,2}.fastq.gz (kraken2 splits on the required # placeholder); single-end → &lt;sample&gt;_unclassified.fastq.gz. Skipped by default to save disk.'
    @params['save_unclassified', "context"] = "Kraken"
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'specify other commandline options for kraken; do not specify any option that is already covered by the dedicated input fields'
    @params['cmdOptions', "context"] = "Kraken"
    
        # trimming options
    # general
    @params['trimAdapter', 'hr'] = true
    @params['trimAdapter'] = true
    @params['trimAdapter', "context"] = "OpenGene/fastp"
    # Fastp
    ## trimming
    @params['trim_front1'] = '0'
    @params['trim_front1','description'] = 'trimming how many bases in front for read1 (and read2), default is 0.'
    @params['trim_front1', "context"] = "OpenGene/fastp"
    @params['trim_tail1'] = '0'
    @params['trim_tail1','description'] = 'trimming how many bases in tail for read1 (and read2), default is 0.'
    @params['trim_tail1', "context"] = "OpenGene/fastp"
    @params['cut_front'] = false
    @params['cut_front','description'] = 'move a sliding window from front (5p) to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.'
    @params['cut_front', "context"] = "OpenGene/fastp"
    @params['cut_front_window_size'] = '4'
    @params['cut_front_window_size', "context"] = "OpenGene/fastp"
    @params['cut_front_mean_quality'] = '20'
    @params['cut_front_mean_quality', "context"] = "OpenGene/fastp"
    @params['cut_tail'] = false
    @params['cut_tail','description'] = 'move a sliding window from tail (3p) to front, drop the bases in the window if its mean quality < threshold, stop otherwise.'
    @params['cut_tail', "context"] = "OpenGene/fastp"
    @params['cut_tail_window_size'] = '4'
    @params['cut_tail_window_size', "context"] = "OpenGene/fastp"
    @params['cut_tail_mean_quality'] = '20'
    @params['cut_tail_mean_quality', "context"] = "OpenGene/fastp"
    @params['cut_right'] = false
    @params['cut_right','description'] = 'move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, and then stop.'
    @params['cut_right', "context"] = "OpenGene/fastp"
    @params['cut_right_window_size'] = '4'
    @params['cut_right_window_size', "context"] = "OpenGene/fastp"
    @params['cut_right_mean_quality'] = '20'
    @params['cut_right_mean_quality', "context"] = "OpenGene/fastp"
    @params['average_qual'] = '0'
    @params['average_qual','description'] = 'if one read average quality score <avg_qual>, then this read/pair is discarded. Default 0 means no requirement'
    @params['average_qual', "context"] = "OpenGene/fastp"
    @params['max_len1'] = '0'
    @params['max_len1','description'] = 'if read1 is longer than max_len1, then trim read1 at its tail to make it as long as max_len1. Default 0 means no limitation. If two reads are present, the same will apply to read2.'
    @params['max_len1', "context"] = "OpenGene/fastp"
    @params['max_len2'] = '0'
    @params['max_len2','description'] = 'if read1 is longer than max_len2, then trim read2 at its tail to make it as long as max_len1. Default 0 means no limitation.'
    @params['max_len2', "context"] = "OpenGene/fastp"
    @params['poly_x_min_len'] = '10'
    @params['poly_x_min_len','description'] = 'the minimum length to detect polyX in the read tail. 10 by default.'
    @params['poly_x_min_len', "context"] = "OpenGene/fastp"
    @params['length_required'] = '18'
    @params['length_required','description'] = 'reads shorter than length_required will be discarded.'
    @params['length_required', "context"] = "OpenGene/fastp"
    @params['cmdOptionsFastp'] = ''
    @params['cmdOptionsFastp', "context"] = "OpenGene/fastp"
    @params['mail'] = ""
    @modules = ["QC/fastp", "Dev/R", "Tools/kraken/2.17.1", "Tools/Krona"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def preprocess
    # 'exclusive' drives the whole fast path: it both reserves the node (via the
    # framework's submit()) and, here, switches the run to a single dataset job so
    # the samples can share one k2 daemon (see initialize). Deriving process_mode
    # from it keeps the two in lock-step — they can never be set inconsistently.
    # Runs before set_output_files/main and after the GUI has applied user params
    # (submit_job.rb sets @params before run), so the derived value is in effect.
    @params['process_mode'] = (@params['exclusive'].to_s == 'true') ? 'DATASET' : 'SAMPLE'
    if @params['paired']
      @required_columns << 'Read2'
    end
  end
  # Build the output-table row for a single sample.
  # Note: per-read .txt.gz and <name>_unclassified.fasta.gz are written to
  # the result dir but intentionally not surfaced here — only the aggregate
  # report and Krona are advertised downstream.
  def output_row_for(name, read1 = nil, read2 = nil)
    # Krona is now rendered with ktImportText (see ezRun app-kraken.R), which
    # writes a single self-contained <name>.html — there is no <name>.html.files
    # dir anymore, so KronaOutDir is intentionally not emitted (copying a missing
    # dir would fail the g-req step).
    out = {'Name'=>name,
     'KronaReport [Link]'=>File.join(@result_dir, "#{name}.html"),
     'KrakenReport [File]'=>File.join(@result_dir, "#{name}.report.txt"),
     'KronaOut [File]'=>File.join(@result_dir, "#{name}.html"),
     'Live Report [Link]'=>"http://fgcz-shiny.uzh.ch/exploreMetaTax?data=#{@result_dir}",
    }
    # Carry raw FASTQ paths forward so downstream apps that need the reads
    # (Bracken -> HUMAnN, etc.) can match this dataset without a manual
    # merge. Use [Link] not [File]: these are pass-through pointers to the
    # upstream FASTQ folder, NOT produced by this Kraken run. [File] would
    # tell sushiApp.rb:613 to g-req copy the FASTQs back into their own
    # source folder — self-referential copy fails with "destination path
    # already exists" and (via set -e) kills the job. sushi_fabric.rb:338
    # normalises the tag away for required-columns matching, so
    # `Read1 [Link]` still satisfies downstream apps declaring `Read1`.
    out['Read1 [Link]'] = read1 if read1
    out['Read2 [Link]'] = read2 if read2
    out.merge(extract_columns(@inherit_tags, sample_name: name))
  end

  # True if this sample should be processed (honours the optional `samples`
  # re-run filter, same rule the framework's sample_mode/dataset_mode use).
  def selected_sample?(name)
    return true if @params['samples'].to_s.empty?
    @params['samples'].split(',').include?(name)
  end

  # One output-table row per selected sample. In DATASET mode @dataset is an
  # Array of input row hashes (set by set_input_dataset / dataset_mode). We
  # re-apply selected_sample? here so this returns the SAME set regardless of
  # whether @dataset has already been filtered (see next_dataset caveat below).
  def dataset_rows
    samples = @dataset.is_a?(Array) ? @dataset : [@dataset]
    samples.map { |raw| Hash[*raw.map { |k, v| [k.to_s.gsub(/\[.+\]/, '').strip, v] }.flatten] }
           .select { |row| selected_sample?(row['Name']) }
           .map { |row| output_row_for(row['Name'], row['Read1'], row['Read2']) }
  end

  # In SAMPLE mode (default) the framework sets @dataset to one tag-stripped row
  # hash and expects one clean output row back — classic behaviour, unchanged.
  #
  # In DATASET mode (exclusive on) one job covers the whole dataset, and the
  # framework calls next_dataset in three places, none of which is the output
  # table (that is built row-per-sample in dataset_mode below):
  #   (a) set_output_files  -> to learn which headers are [File];
  #   (b) job_footer        -> to know which files to g-req copy scratch->gstore;
  #   (c) run_RApp          -> serialised into the R `output`, which EzApp$run
  #                            turns into an EzDataset (so it needs a 'Name').
  # A single row hash can name only ONE sample's files, which would drop samples
  # 2..N from the copy (b). So we return a *copy manifest*: every selected sample's
  # [File] outputs under unique synthetic keys, plus a 'Name' so (c) can build a
  # valid EzDataset. The worker never reads `output`, so Name's value is cosmetic;
  # these synthetic keys never reach dataset.tsv.
  def next_dataset
    if @params['process_mode'] == 'SAMPLE'
      return output_row_for(@dataset['Name'], @dataset['Read1'], @dataset['Read2'])
    end
    rows = dataset_rows
    manifest = {}
    manifest['Name'] = rows.first['Name'] unless rows.empty?
    rows.each do |row|
      name = row['Name']
      row.each do |header, value|
        next unless header.to_s =~ /\[File\]/
        base = header.to_s.sub(/\s*\[.+\]\s*/, '').strip
        manifest["#{name} #{base} [File]"] = value
      end
    end
    manifest
  end

  # Keep the historical output table (one clean row per sample) even though the
  # framework's dataset_mode collapses next_dataset to a single row. Let the
  # framework build @dataset + the single job script + the file-copy manifest,
  # then replace @result_dataset with the per-sample rows.
  def dataset_mode
    super
    @result_dataset = dataset_rows
  end

  def commands
    run_RApp("EzAppKraken")
  end

  # Scan /srv/GT/databases/kraken2/ for valid Kraken2 DB folders. A folder is
  # a valid DB if it contains 'hash.k2d'. Returns an Array of [label, value]
  # pairs so the SUSHI multi-select widget (which only renders for Arrays,
  # not Hashes) shows on-disk size while submitting just the folder name.
  # Standard is promoted first when present; the rest are alphabetical.
  # Returns [] if the root doesn't exist or contains no valid DBs.
  def kraken_db_choices
    root = '/srv/GT/databases/kraken2'
    return [] unless Dir.exist?(root)

    dbs = Dir.entries(root)
             .reject { |e| e.start_with?('.') }
             .select { |e| File.exist?(File.join(root, e, 'hash.k2d')) }
             .sort

    preferred = 'Standard'
    dbs = [preferred] + (dbs - [preferred]) if dbs.include?(preferred)

    dbs.map do |name|
      bytes = `du -sb #{File.join(root, name)} 2>/dev/null`.to_i
      gb = (bytes.to_f / (1024**3))
      label = bytes > 0 ? "#{name} (#{format('%.1f', gb)} GB)" : name
      [label, name]
    end
  end
end

if __FILE__ == $0
end
