#!/usr/bin/env ruby
# encoding: utf-8
#
# LitmusProcessApp — SUSHI frontend for the LITMUS Processing workbench step.
#
# Wraps scripts/R/processing_engine.R (9 engines: MZmine 4, MS-DIAL 5, XCMS,
# OpenMS, MASSter, MassCube, Skyline, eMZed, El-MAVEN + LipidOracle/Goslin
# annotation). Takes a per-sample dataset.tsv (one row per .raw / .mzML file),
# produces a single feature intensity matrix + metadata CSVs + HTML report.
#
# Design principles adopted from reviewer feedback (2026-04-24):
#   (1) Sample grouping is METADATA carried in dataset.tsv (not parsed from
#       filenames). Columns: "Grouping [Factor]", "Sample Type [Factor]",
#       "Polarity [Factor]". The R script reads these verbatim.
#   (2) Input contract follows the SUSHI "one file per sample" convention
#       (cf. NanoPlotApp, BismarkApp). Output contract consolidates all
#       per-sample information into a single feature_table.csv
#       (samples × features), matching the LITMUS internal data model.
#   (3) Closer to ezRun style than the 2026-04-23 LitmusNormalize PoC:
#       numbered params, cores/ram/scratch with 'slurm' context, inherit
#       Factor / B-Fabric tags, grouping-var + sample-type gates exposed
#       as dropdowns so the UI reflects the experimental design.
#
# Branch: litmus-sushi-app
# PoC instance: /srv/sushi/masa_test_sushi_20260416/

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class LitmusProcessApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'LitmusProcess'
    @analysis_category = 'Lipidomics'
    @description = <<-EOS
LITMUS Processing — untargeted LC-MS feature detection across a cohort of .raw
or .mzML files. Runs one of six feature-detection engines (XCMS, MZmine 4,
MS-DIAL 5, OpenMS, MASSter, MassCube) over all selected samples, produces a
single <b>samples × features</b> intensity matrix plus feature &amp; sample
metadata and an HTML report.<br/>
Input: SUSHI dataset.tsv with one row per acquisition file. Sample grouping
(<i>Grouping</i>, <i>Sample Type</i>, <i>Polarity</i>) must be present as
metadata columns — the app does <b>not</b> parse filenames for classification.<br/>
Backend: <code>scripts/R/processing_engine.R</code>
(<a href='https://github.com/masaomi/LITMUS'>LITMUS core</a>).
EOS

    @params['process_mode'] = 'DATASET'

    # --- SUSHI input contract ---------------------------------------------
    # Required columns (SUSHI strips the [File] / [Factor] type suffix when
    # matching) — these are the ezRun-style canonical column names.
    @required_columns = ['Name', 'Thermo Scientific RAW', 'Polarity', 'Sample Type']
    @required_params  = ['engine', 'polarity_filter', 'sample_types']

    # --- SLURM resources --------------------------------------------------
    # XCMS with 50-150 Orbitrap files needs ~16-32 GB RAM and can soak up
    # all cores on peak detection. Leave generous defaults; users override
    # for small cohorts.
    @params['cores']   = '16'
    @params['cores',   'context'] = 'slurm'
    @params['ram']     = '32'
    @params['ram',     'context'] = 'slurm'
    @params['scratch'] = '100'
    @params['scratch', 'context'] = 'slurm'

    # --- Engine selection (numbered like DIANN for stable UI ordering) ----
    @params['01_engine'] = ['xcms', 'mzmine4', 'msdial5', 'openms', 'masster', 'masscube']
    @params['01_engine', 'description'] = 'Feature-detection engine. XCMS is the default — pure Bioconductor R, no external binary needed.'

    # Convenience alias (needed in @required_params)
    @params['engine'] = 'xcms'
    @params['engine', 'description'] = 'Alias of 01_engine (kept for backwards compatibility with pipelines that reference the unversioned key).'

    # --- Sample gating (polarity, sample-type) ----------------------------
    # These filter which rows of the dataset.tsv actually enter the pipeline.
    # Per reviewer feedback: filtering by metadata, never by filename parsing.
    @params['02_polarity_filter'] = ['pos', 'neg', 'both']
    @params['02_polarity_filter', 'description'] = 'Which polarity to process. Separate jobs recommended for pos / neg.'
    @params['polarity_filter']    = 'pos'

    @params['03_sample_types_include'] = 'sample,qc'
    @params['03_sample_types_include', 'description'] = 'Comma-separated list of Sample Type values to INCLUDE. blanks/standards are typically excluded from untargeted feature detection but kept in the sample_meta.csv for traceability.'
    @params['sample_types']       = 'sample,qc'

    # --- Untargeted peak detection defaults (Orbitrap) --------------------
    # All parameters in ppm or minutes per CLAUDE.md instrument rules.
    @params['11_mz_tolerance_ppm']  = '5'
    @params['11_mz_tolerance_ppm',  'description'] = 'MS1 mass tolerance (ppm). 5 ppm is the LITMUS Orbitrap default (safety margin above 1-3 ppm instrument spec).'

    @params['12_rt_peakwidth_sec']  = '5,40'
    @params['12_rt_peakwidth_sec',  'description'] = 'CentWave peak width range (seconds). 5-40 s covers sub-2 µm UHPLC.'

    @params['13_snthresh']          = '3'
    @params['13_snthresh',          'description'] = 'CentWave signal-to-noise threshold.'

    @params['14_noise_level']       = '1000'
    @params['14_noise_level',       'description'] = 'Intensity floor for peak detection.'

    @params['15_min_frac']          = '0.5'
    @params['15_min_frac',          'description'] = 'Minimum fraction of samples a feature must appear in (post-grouping).'

    # --- Lipid annotation (post-processing) -------------------------------
    @params['21_lipid_annotation']  = ['none', 'goslin', 'lipidoracle']
    @params['21_lipid_annotation',  'description'] = 'Optional lipid-name standardization / annotation applied to feature_meta.csv after detection. Goslin is pure R; LipidOracle requires Docker.'

    # --- Demo / offline mode ----------------------------------------------
    # When the dataset.tsv points to .raw files that don't exist on the
    # compute node (PoC scenario), the R script generates a synthetic
    # feature matrix so the SUSHI wiring (epilogue, report render, gstore
    # copy) can be end-to-end validated without real mass-spec input.
    @params['91_demo_mode']         = ['auto', 'always', 'never']
    @params['91_demo_mode',         'description'] = 'auto: fall back to synthetic data if fewer than 2 real files survive polarity + sample-type filtering (engines need >=2 files). always: synthetic regardless of input. never: error out if any expected file is missing.'

    # --- Paths (shared analysis + LITMUS library) -------------------------
    @params['99_litmus_path']       = '/srv/GT/analysis/masaomi/2026/FGCZ/litmus_sushiapp_process_20260424'
    @params['99_litmus_path',       'description'] = 'Staged LITMUS scripts directory (must contain processing_engine.R + litmus_process.R + Rmd template).'

    @params['name'] = 'LitmusProcess_Result'
    @params['mail'] = ''

    @modules       = ['Dev/R']
    @inherit_tags  = ['Factor', 'B-Fabric']
  end

  # SUSHI next_dataset: the parent directory is the single gstore drop-off;
  # individual files are declared as [Link] so they render in the UI but
  # don't each trigger a separate g-req copy (cf. LitmusNormalize PoC
  # gotcha #4, 2026-04-23).
  def next_dataset
    rel_dir = File.join(@result_dir, @params['name'])
    {
      'Name'                     => @params['name'],
      'Result [File]'            => rel_dir,
      'Report [Link]'            => File.join(rel_dir, '00index.html'),
      'Feature Table [Link]'     => File.join(rel_dir, 'feature_table.csv'),
      'Feature Meta [Link]'      => File.join(rel_dir, 'feature_meta.csv'),
      'Sample Meta [Link]'       => File.join(rel_dir, 'sample_meta.csv'),
      'Processing Log [Link]'    => File.join(rel_dir, 'processing.log'),
      'engine'                   => @params['01_engine'],
      'polarity_filter'          => @params['02_polarity_filter']
    }.merge(extract_columns(@inherit_tags))
  end

  def commands
    script_path  = File.join(@params['99_litmus_path'], 'scripts', 'R', 'litmus_process.R')
    out_dir      = @params['name']

    # @input_dataset_tsv_path is SUSHI's staged copy of the submission
    # dataset.tsv on the compute node (see global_variables.rb line 328).
    dataset_tsv  = @input_dataset_tsv_path

    # Resolve paths inside dataset.tsv at R-time — the [File] column holds
    # gstore-relative paths (e.g. "p35611/Metabolomics/EXPLORIS_4/..."),
    # R prepends the gstore base.

    <<~CMD.strip
      mkdir -p '#{out_dir}'
      Rscript --vanilla '#{script_path}' \\
        --dataset '#{dataset_tsv}' \\
        --output  '#{out_dir}' \\
        --engine  '#{@params['01_engine']}' \\
        --polarity '#{@params['02_polarity_filter']}' \\
        --sample_types '#{@params['03_sample_types_include']}' \\
        --mz_tol_ppm '#{@params['11_mz_tolerance_ppm']}' \\
        --rt_peakwidth_sec '#{@params['12_rt_peakwidth_sec']}' \\
        --snthresh '#{@params['13_snthresh']}' \\
        --noise '#{@params['14_noise_level']}' \\
        --min_frac '#{@params['15_min_frac']}' \\
        --annotation '#{@params['21_lipid_annotation']}' \\
        --demo_mode '#{@params['91_demo_mode']}' \\
        --litmus_path '#{@params['99_litmus_path']}' \\
        --gstore_base /srv/gstore/projects
    CMD
  end
end

if __FILE__ == $0
  # CLI invocation for PoC testing:
  #   cd /srv/sushi/masa_test_sushi_20260416/master
  #   bundle exec sushi_fabric --class LitmusProcessApp \
  #     --dataset /srv/sushi/masa_test_sushi_20260416/test_data/litmus_process_dataset.tsv \
  #     --project 35611 --run
end
