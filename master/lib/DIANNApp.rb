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
3-step snakemake (lib search -> quant refinement -> final quant) with
ThermoRawFileParser conversion and prolfquapp QC. Outputs a Result_WU&lt;id&gt;.zip
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
    @required_params = ['name', 'paramsTemplate']

    # ---- Compute + identity ----
    @params['cores']   = '64'   ; @params['cores',   'context']  = 'slurm'
    @params['ram']     = '32'   ; @params['ram',     'context']  = 'slurm'
    @params['scratch'] = '20'   ; @params['scratch', 'context']  = 'slurm'
    @params['node']    = ['fgcz-c-050']
    @params['node', 'context'] = 'slurm'
    @params['name']    = 'DIANN_v23'
    @params['mail']    = ''

    # ---- Template chooser ----
    # paramsTemplate sets the BASE for the params block; per-key fields below
    # override it, and customParamsYml replaces everything.
    @params['paramsTemplate']  = ['default-DIA', 'default-DDA']
    @params['customParamsYml'] = ''       # path; full override when set

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

    # ---- DIA / DDA mode ----
    @params['is_dda']                         = ['false', 'true']
    @params['scan_window']                    = ['AUTO', '7', '11', '15']

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
    @params['mass_acc_ms1']                   = ['AUTO', '5', '10', '15', '20']
    @params['mass_acc_ms2']                   = ['AUTO', '5', '10', '15', '20']

    # ---- Scoring + protein inference ----
    @params['scoring_qvalue']                 = ['0.01', '0.001', '0.05']
    @params['protein_pg_level']               = ['protein_names_1', 'genes_0', 'isoforms_2']
    @params['protein_relaxed_prot_inf']       = ['false', 'true']

    # ---- Quantification ----
    @params['quantification_reanalyse']       = ['true', 'false']
    @params['quantification_no_norm']         = ['false', 'true']

    # ---- Freestyle escape hatch ----
    @params['freestyle']                      = 'None'  # raw DIA-NN flags appended verbatim

    # ---- Conversion + verbosity ----
    @params['raw_converter']                  = ['thermoraw', 'msconvert', 'msconvert-demultiplex']
    @params['verbose']                        = ['1', '0', '2', '3']

    @modules = ['Dev/R']  # snakemake env activated via @conda_env below
    @inherit_columns = []
  end

  def next_dataset
    # [File] basenames MUST match the dirs run-diann produces in the job cwd —
    # SUSHI's g-req copies `File.basename(value)` from cwd to <gstore>/<dirname>.
    # run-diann writes `qc_result/` and `out-DIANN_quantC/` directly, so no
    # rename step is needed in the app.
    {
      'Name'                      => @params['name'],
      'Protein Abundances [Link]' => File.join(@result_dir, 'qc_result/proteinAbundances.html'),
      'Sample Sizes [Link]'       => File.join(@result_dir, 'qc_result/QC_sampleSizeEstimation.html'),
      'qc_result [File]'          => File.join(@result_dir, 'qc_result'),
      'DIANN Quant [File]'        => File.join(@result_dir, 'out-DIANN_quantC')
    }.merge(extract_columns(colnames: @inherit_columns))
  end

  def commands
    # Python entry point: from ezpyz_diann.app import EzAppDiann
    # (run_PyApp resolves Ruby app name 'DIANN' to module ezpyz_diann + class EzAppDiann)
    run_PyApp('DIANN', conda_env: 'gi_snakemake8.20.5')
  end
end

if __FILE__ == $0

end
