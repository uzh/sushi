#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class DIANNApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'DIANN'
    @analysis_category = 'Proteomics'
    @description =<<-EOS
DIA / DDA proteomics quantification with DIA-NN 2.3 + prolfqua QC.<br/>
Pipeline: <a href="https://github.com/wolski/diann-runner">diann-runner</a>
two-step snakemake by default (lib search -> quant refinement), with optional
Step C final quantification, ThermoRawFileParser conversion, and prolfquapp QC.
Outputs a Result_WU&lt;id&gt;.zip
plus interactive proteinAbundances.html / QC_sampleSizeEstimation.html.<br/>
Container runtime is apptainer (rootless), so the job runs as trxcopy on a
genomics compute node.
    EOS
    @params['process_mode'] = 'DATASET'
    # 'Thermo RAW' scopes this app to proteomics datasets (B-Fabric standard
    # column for Thermo raw mass spec files). EzAppDiann ALSO accepts the
    # legacy 'RAW' column as a fallback for CLI submissions or non-B-Fabric
    # datasets, but the SUSHI UI gates on the standard column.
    @required_columns = ['Name', 'Thermo RAW']
    @required_params = ['name', 'diann_version']

    # ---- Compute + identity ----
    @params['cores']   = ['8','2']   ; @params['cores',   'context']  = 'slurm'
    @params['ram']     = '32'   ; @params['ram',     'context']  = 'slurm'
    @params['scratch'] = '20'   ; @params['scratch', 'context']  = 'slurm'
    @params['node']    = ['fgcz-c-050']
    @params['node', 'context'] = 'slurm'
    @params['name']    = 'DIANN_v23'
    @params['mail']    = ''

    # The GUI defaults below ARE the known-good DIA preset (provenance: FGCZ
    # workunit WU340602 / p40993, a DIA-NN 2.3.2 DIA run on Orbitrap Exploris
    # data). run-diann maps every readable key directly onto its internal params
    # (SUSHI_TO_DRUNNER) — there is no B-Fabric-keyed template indirection, so the
    # values here are the single source. For DDA, set is_dda + the mass-acc /
    # charge fields accordingly.

    # ---- DIA-NN version ----
    # Selects the DIA-NN container image: run-diann maps this value to
    # diann_images[<version>]. First entry is the GUI default.
    @params['diann_version'] = ['2.5.1', '2.5.0', '2.3.2']

    # ---- FASTA databases ----
    # order_fasta: checkbox for now — when checked, use the order-specific FASTA.
    @params['order_fasta'] = false
    # fasta_databases: multi-select populated by globbing the FASTA library on
    # the SUSHI host (read at form-render time, so FASTA_DB_DIR must be mounted
    # there). Options are full paths; selected entries come back comma-joined and
    # the app passes one --fasta per entry to run-diann.
    fasta_db_dir = '/srv/GT/databases/proteomics/fasta'  # TODO: set to the real FASTA library on the SUSHI host
    @params['fasta_databases'] = Dir.glob(File.join(fasta_db_dir, '*.{fasta,fa}')).sort
    @params['fasta_databases', 'multi_selection'] = true

    # ---- Workflow ----
    @params['workflow_mode'] = ['two_step', 'one_step']
    @params['enable_step_c'] = ['false', 'true']

    # ---- DIA / DDA mode ----
    @params['is_dda']                         = ['false', 'true']
    @params['scan_window']                    = ['AUTO', '7', '11', '15']
    # "Unrelated runs" = --individual-mass-acc --individual-windows (per-run calibration)
    @params['unrelated_runs']                 = ['false', 'true']

    # ---- Modifications ----
    @params['mods_variable']                  = '--var-mods 1 --var-mod UniMod:35,15.994915,M'
    @params['mods_no_peptidoforms']           = ['false', 'true']
    @params['mods_unimod4']                   = ['true', 'false']
    @params['mods_met_excision']              = ['true', 'false']

    # ---- Peptide constraints ----
    @params['peptide_min_length']             = '6'
    @params['peptide_max_length']             = '30'
    @params['peptide_precursor_charge_min']   = '1'
    @params['peptide_precursor_charge_max']   = '3'
    @params['peptide_precursor_mz_min']       = '350'
    @params['peptide_precursor_mz_max']       = '1500'
    @params['peptide_fragment_mz_min']        = '200'
    @params['peptide_fragment_mz_max']        = '1800'

    # ---- Digestion ----
    @params['digestion_cut']                  = ['K*,R*', 'K*', 'R*']
    @params['digestion_missed_cleavages']     = ['1', '0', '2', '3']

    # ---- Mass accuracy ----
    # See TODO_mass_accuracy: AUTO default; 4/7 added for high-res Orbitrap/Astral.
    @params['mass_acc_ms1']                   = ['AUTO', '4', '5', '7', '10', '15', '20']
    @params['mass_acc_ms2']                   = ['AUTO', '4', '5', '7', '10', '15', '20']

    # ---- Scoring + protein inference ----
    @params['scoring_qvalue']                 = ['0.01', '0.001', '0.05']
    @params['protein_pg_level']               = ['2_genes', '1_protein_names', '0_isoform_IDs']

    # ---- Quantification ----
    @params['quantification_reanalyse']       = ['true', 'false']
    @params['quantification_no_norm']         = ['false', 'true']
    # Export fragment-level per-run quantities (--export-quant); off by default
    # as it substantially enlarges the report.
    @params['quantification_export_quant']    = ['false', 'true']

    # ---- Freestyle escape hatch ----
    @params['freestyle']                      = 'None'  # raw DIA-NN flags appended verbatim

    # ---- Conversion + verbosity ----
    @params['raw_converter']                  = ['native','thermoraw', 'msconvert', 'msconvert-demultiplex']
    @params['verbose']                        = ['1', '0', '2', '3']

    @modules = ['Dev/R', 'Dev/pixi']  # snakemake env activated via @conda_env below; pixi runs the DIANN app
    @inherit_columns = ['Order Id', 'Thermo RAW']
  end

  def next_dataset
    # [File] basenames MUST match the dirs run-diann produces in the job cwd —
    # SUSHI's g-req copies `File.basename(value)` from cwd to <gstore>/<dirname>.
    # run-diann writes `qc_result/` and `out-DIANN_quantB/quantC` directly, so no
    # rename step is needed in the app.
    quant_dir = @params['enable_step_c'].to_s == 'true' ? 'out-DIANN_quantC' : 'out-DIANN_quantB'
    {
      'Name'                      => @params['name'],
      'Protein Abundances [Link]' => File.join(@result_dir, 'qc_result/proteinAbundances.html'),
      'Sample Sizes [Link]'       => File.join(@result_dir, 'qc_result/QC_sampleSizeEstimation.html'),
      'qc_result [File]'          => File.join(@result_dir, 'qc_result'),
      'DIANN Quant [File]'        => File.join(@result_dir, quant_dir)
    }.merge(extract_columns(colnames: @inherit_columns))
  end

  def commands
    # Python entry point: from ezpyz_diann.app import EzAppDiann
    # (run_PyApp resolves Ruby app name 'DIANN' to module ezpyz_diann + class EzAppDiann)
    run_PyApp('DIANN', pixi_enabled: true)  # Name must match [name] in 'ezpyz_[name]' format
  end
end

if __FILE__ == $0

end
