#!/usr/bin/env ruby
# encoding: utf-8
#
# LitmusNormalizeApp — SUSHI frontend for LITMUS lipidomics normalization.
#
# Wraps normalize_lipidomics.R (LITMUS core library, 120+ functions, matrix-in/matrix-out).
# Accepts a feature intensity matrix (CSV, samples-in-rows), runs a configurable
# normalize → transform → scale pipeline, outputs normalized CSV + QC metrics + HTML report.
#
# PoC branch: litmus-sushi-app (fgcz-h-083 test instance)

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class LitmusNormalizeApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'LitmusNormalize'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Lipidomics'
    @description = <<-EOS
LITMUS Lipidomics Normalization — applies sample-level normalization (PQN / TIC / Median / Quantile / IS),
optional log transform (log2 / log10), and feature-level scaling (Pareto / Auto) to a feature intensity matrix.<br/>
Input: CSV with samples in rows, features in columns (first column = sample name, optional Group column).<br/>
Output: Normalized matrix CSV + QC metrics CSV + HTML report.<br/>
Core library: <a href='https://github.com/masaomi/LITMUS'>LITMUS normalize_lipidomics.R</a>
EOS
    @required_columns = ['Name', 'Intensity']
    @required_params  = ['normalization_method', 'transform', 'scaling']

    # Computational resources — small for lipidomics matrices (~100 samples × ~1000 features)
    @params['cores']   = '2'
    @params['cores',   'context'] = 'slurm'
    @params['ram']     = '8'
    @params['ram',     'context'] = 'slurm'
    @params['scratch'] = '10'
    @params['scratch', 'context'] = 'slurm'

    # LITMUS-specific parameters
    @params['normalization_method'] = ['PQN', 'TIC', 'Median', 'Quantile', 'IS', 'none']
    @params['normalization_method', 'description'] = 'Sample-level normalization method'

    @params['transform'] = ['log2', 'log10', 'none']
    @params['transform', 'description'] = 'Mathematical transformation applied after normalization'

    @params['scaling'] = ['pareto', 'auto', 'none']
    @params['scaling', 'description'] = 'Feature-level scaling applied after transform'

    @params['groups_file'] = ''
    @params['groups_file', 'description'] = '(Optional) absolute path to group metadata CSV with SampleID,Group columns. Leave blank to derive groups from the Group column of the intensity CSV.'

    @params['litmus_path'] = '/srv/GT/analysis/masaomi/2026/FGCZ/litmus_sushiapp_test_20260423'
    @params['litmus_path', 'description'] = 'Path to LITMUS repository root (provides normalize_lipidomics.R)'

    @params['name'] = 'LitmusNormalize_Result'
    @params['mail'] = ''

    @modules = ['Dev/R']
    @inherit_tags = ['Factor', 'B-Fabric']
  end

  def next_dataset
    # SUSHI convention (ezRun-style):
    # - Result [File] points to the whole OUTPUT DIRECTORY.
    #   g-req will copy that directory (created under scratch) to gstore.
    # - Individual output files are declared as [Link] so their URLs render
    #   in the Web UI after the directory is served from gstore, but they
    #   don't trigger a separate g-req copy (which would fail — the files
    #   live inside the directory, not at scratch root).
    rel_dir = File.join(@result_dir, @params['name'])
    {
      'Name'                 => @params['name'],
      'Result [File]'        => rel_dir,
      'Report [Link]'        => File.join(rel_dir, '00index.html'),
      'Normalized [Link]'    => File.join(rel_dir, 'normalized_matrix.csv'),
      'QC Metrics [Link]'    => File.join(rel_dir, 'qc_metrics.csv'),
      'normalization_method' => @params['normalization_method'],
      'transform'            => @params['transform'],
      'scaling'              => @params['scaling']
    }.merge(extract_columns(@inherit_tags))
  end

  def commands
    # R script and normalize_lipidomics.R both live under <litmus_path>/scripts/R/
    script_path = File.join(@params['litmus_path'], 'scripts', 'R', 'litmus_normalize.R')

    # Resolve a dataset [File] value that may be either:
    #   (a) absolute path (starts with /)  — use verbatim
    #   (b) gstore-relative path (e.g. "p35611/.../file.csv") — prepend gstore base
    gstore_base = '/srv/gstore/projects'
    resolve = ->(p) {
      s = p.to_s.strip
      next nil if s.empty?
      s.start_with?('/') ? s : File.join(gstore_base, s)
    }

    input_file  = resolve.call(@dataset[0]['Intensity [File]'])
    groups_file = resolve.call(@dataset[0]['Groups [File]']) \
                  || (resolve.call(@params['groups_file']) rescue nil)
    # Output goes into a subdir of the scratch working dir. SUSHI's epilogue
    # does `g-req -w copy <subdir> <gstore_dest_parent>` (see next_dataset's
    # Result [File]) — so we must NOT write to /srv/gstore/ directly (read-only
    # on compute nodes) and NOT prepend p35611/<job>/ here (that gets created
    # implicitly by SUSHI's dest_dir derivation).
    out_dir     = @params['name']
    groups_arg  = groups_file ? " \\\n        --groups '#{groups_file}'" : ''

    <<~CMD.strip
      mkdir -p '#{out_dir}'
      Rscript --vanilla '#{script_path}' \\
        --input '#{input_file}' \\
        --output '#{out_dir}' \\
        --method '#{@params['normalization_method']}' \\
        --transform '#{@params['transform']}' \\
        --scaling '#{@params['scaling']}' \\
        --litmus_path '#{@params['litmus_path']}'#{groups_arg}
    CMD
  end
end

if __FILE__ == $0
  # Enable CLI invocation via sushi_fabric for PoC testing:
  #   cd /srv/sushi/<test>/master
  #   bundle exec sushi_fabric --class LitmusNormalizeApp \
  #     --dataset /srv/sushi/<test>/dataset.tsv \
  #     --project 1535 --run
end
