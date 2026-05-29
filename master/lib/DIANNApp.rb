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
    # paramsTemplate sets the BASE for the params block. Per-key fields
    # below override it; customParamsYml replaces everything.
    @params['paramsTemplate']  = ['default-DIA', 'default-DDA']
    @params['customParamsYml'] = ''       # path; full override when set
    @params['order_fasta']     = ''       # text path; staged as input/order.fasta

    # ---- Workflow ----
    @params['02_workflow_mode']                    = ['two_step', 'one_step']
    @params['03_fasta_database_path']              = 'NONE'           # text path
    @params['03_fasta_use_custom']                 = ['false', 'true']
    @params['03b_additional_fasta_database_path']  = 'NONE'           # text path

    # ---- DIA / DDA mode ----
    @params['05_diann_is_dda']                     = ['false', 'true']
    @params['05b_diann_scan_window']               = ['AUTO', '7', '11', '15']

    # ---- Modifications ----
    @params['06a_diann_mods_variable']             = '--var-mods 1 --var-mod UniMod:35,15.994915,M'
    @params['06b_diann_mods_no_peptidoforms']      = ['false', 'true']
    @params['06c_diann_mods_unimod4']              = ['true', 'false']
    @params['06d_diann_mods_met_excision']         = ['true', 'false']

    # ---- Peptide constraints ----
    @params['07_diann_peptide_min_length']           = '6'
    @params['07_diann_peptide_max_length']           = '30'
    @params['07_diann_peptide_precursor_charge_min'] = '1'
    @params['07_diann_peptide_precursor_charge_max'] = '3'
    @params['07_diann_peptide_precursor_mz_min']     = '350'
    @params['07_diann_peptide_precursor_mz_max']     = '1500'
    @params['07_diann_peptide_fragment_mz_min']      = '200'
    @params['07_diann_peptide_fragment_mz_max']      = '1800'

    # ---- Digestion ----
    @params['08_diann_digestion_cut']              = ['K*,R*', 'K*', 'R*']
    @params['08_diann_digestion_missed_cleavages'] = ['1', '0', '2', '3']

    # ---- Mass accuracy ----
    @params['09_diann_mass_acc_ms1']               = ['AUTO', '5', '10', '15', '20']
    @params['09_diann_mass_acc_ms2']               = ['AUTO', '5', '10', '15', '20']

    # ---- Scoring + protein inference ----
    @params['10_diann_scoring_qvalue']             = ['0.01', '0.001', '0.05']
    @params['11a_diann_protein_pg_level']          = ['protein_names_1', 'genes_0', 'isoforms_2']
    @params['11b_diann_protein_relaxed_prot_inf']  = ['false', 'true']

    # ---- Quantification ----
    @params['12a_diann_quantification_reanalyse']  = ['true', 'false']
    @params['12b_diann_quantification_no_norm']    = ['false', 'true']

    # ---- Freestyle escape hatch ----
    @params['13_diann_freestyle']                  = 'None'  # raw DIA-NN flags appended verbatim

    # ---- Conversion + verbosity ----
    @params['97_raw_converter']  = ['thermoraw', 'msconvert', 'msconvert-demultiplex']
    @params['99_other_verbose']  = ['1', '0', '2', '3']

    @modules = ['Dev/R']  # snakemake env activated via @conda_env below
    @inherit_columns = []
  end

  def next_dataset
    {
      'Name'                      => @params['name'],
      'Protein Abundances [Link]' => File.join(@result_dir, 'qc_result/proteinAbundances.html'),
      'Sample Sizes [Link]'       => File.join(@result_dir, 'qc_result/QC_sampleSizeEstimation.html'),
      'qc_result [File]'          => File.join(@result_dir, 'qc_result'),
      'DIANN Quant [File]'        => File.join(@result_dir, 'DIANN_quantC')
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
