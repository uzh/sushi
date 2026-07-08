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

    # ---- Pipeline / run setup ----
    # diann_version selects the DIA-NN container image: run-diann maps this value
    # to diann_images[<version>]. First entry is the GUI default.
    @params['pipeline_diann_version']         = ['2.5.1', '2.5.0', '2.3.2']
    @params['pipeline_workflow_mode']         = ['two_step', 'single_step']
    @params['pipeline_is_dda']                = ['false', 'true']
    @params['pipeline_raw_converter']         = ['native', 'thermoraw', 'msconvert', 'msconvert-demultiplex']
    @params['enable_step_c']                  = ['false', 'true']

    # ---- Input: FASTA databases ----
    # input_fasta_use_custom: checkbox for now — when checked, use the order-specific FASTA.
    @params['input_fasta_use_custom'] = false
    # input_fasta_databases: multi-select populated by globbing the FASTA library on
    # the SUSHI host (read at form-render time, so FASTA_DB_DIR must be mounted
    # there). Options are full paths; selected entries come back comma-joined and
    # the app passes one --fasta per entry to run-diann.
    fasta_db_dir = '/srv/GT/databases/proteomics/fasta'  # TODO: set to the real FASTA library on the SUSHI host
    @params['input_fasta_databases'] = Dir.glob(File.join(fasta_db_dir, '*.{fasta,fa}')).sort
    @params['input_fasta_databases', 'multi_selection'] = true

    # ---- Library generation: digestion ----
    # 'no digestion' -> empty --cut: a pre-digested / peptide-list FASTA (e.g. the
    # entrapment FASTA) is taken verbatim. Enum aligned with the B-Fabric
    # executable (source of truth): K*,R* / K*,R*,!*P / no digestion.
    @params['lib_digestion_cut']              = ['K*,R*', 'K*,R*,!*P', 'no digestion']
    @params['lib_digestion_missed_cleavages'] = ['1', '0', '2', '3']

    # ---- Library generation: peptide / precursor / fragment ranges ----
    @params['lib_peptide_min_length']         = '6'
    @params['lib_peptide_max_length']         = '30'
    @params['lib_precursor_charge_min']       = '1'
    @params['lib_precursor_charge_max']       = '3'
    @params['lib_precursor_mz_min']           = '350'
    @params['lib_precursor_mz_max']           = '1500'
    @params['lib_fragment_mz_min']            = '200'
    @params['lib_fragment_mz_max']            = '1800'

    # ---- Library generation: modifications ----
    @params['lib_mods_variable']              = '--var-mods 1 --var-mod UniMod:35,15.994915,M'
    @params['lib_mods_no_peptidoforms']       = ['false', 'true']
    @params['lib_mods_unimod4']               = ['true', 'false']
    @params['lib_mods_met_excision']          = ['true', 'false']

    # ---- Search & scoring: mass accuracy ----
    # See TODO_mass_accuracy: AUTO default; 4/7 added for high-res Orbitrap/Astral.
    @params['search_mass_acc_ms1']            = ['AUTO', '4', '5', '7', '10', '15', '20']
    @params['search_mass_acc_ms2']            = ['AUTO', '4', '5', '7', '10', '15', '20']
    # "Unrelated runs" = --individual-mass-acc --individual-windows (per-run calibration)
    @params['search_mass_acc_unrelated_runs'] = ['false', 'true']

    # ---- Search & scoring: FDR + protein inference ----
    @params['search_scoring_qvalue']          = ['0.01', '0.001', '0.05']
    @params['search_protein_pg_level']        = ['2_genes', '1_protein_names', '0_isoform_IDs']
    @params['search_protein_ids_to_names']    = ['false', 'true']

    # ---- Quantification ----
    @params['quant_scan_window']              = ['AUTO', '7', '11', '15']
    @params['quant_reanalyse']                = ['true', 'false']
    @params['quant_no_norm']                  = ['false', 'true']

    # ---- Output ----
    # output_fragment_quant: export fragment-level per-run quantities (--export-quant);
    # off by default as it substantially enlarges the report.
    @params['output_fragment_quant']          = ['false', 'true']
    @params['output_include_libs']            = ['false', 'true']
    @params['output_pmultiqc']                = ['true', 'false']

    # ---- Advanced ----
    @params['advanced_freestyle']             = 'None'  # raw DIA-NN flags appended verbatim
    @params['advanced_verbose']               = ['1', '0', '2', '3']

    @modules = ['Dev/R', 'Dev/pixi']  # snakemake env activated via @conda_env below; pixi runs the DIANN app
    @inherit_columns = ['Order Id', 'Thermo RAW']
    # Carry the sample [Factor] columns forward. DIA-NN produces a single
    # aggregated quant result (one child row), which cannot hold per-sample
    # group assignments. The DEA app needs them, so we ALSO emit one grandchild
    # row per sample (see grandchild_datasets) carrying that sample's [Factor]
    # columns plus the shared quant pointer.
    @inherit_tags = ['Factor']
  end

  def next_dataset
    # [File] basenames MUST match the dirs run-diann produces in the job cwd —
    # SUSHI's g-req copies `File.basename(value)` from cwd to <gstore>/<dirname>.
    # run-diann writes `qc_result/` and `out-DIANN_quantB/quantC` directly, so no
    # rename step is needed in the app.
    quant_dir = @params['enable_step_c'].to_s == 'true' ? 'out-DIANN_quantC' : 'out-DIANN_quantB'
    # extract_columns is `if tags ... elsif colnames`, so the two kinds of
    # inheritance must be merged in separate calls.
    {
      'Name'                      => @params['name'],
      'Protein Abundances [Link]' => File.join(@result_dir, 'qc_result/proteinAbundances.html'),
      'Sample Sizes [Link]'       => File.join(@result_dir, 'qc_result/QC_sampleSizeEstimation.html'),
      'qc_result [File]'          => File.join(@result_dir, 'qc_result'),
      'DIANN Quant [File]'        => File.join(@result_dir, quant_dir)
    }.merge(extract_columns(colnames: @inherit_columns))
     .merge(extract_columns(tags: @inherit_tags))
  end

  # Per-sample grandchild datasets so a downstream DEA app has per-sample group
  # assignments (the aggregated child row cannot). One row per input raw file,
  # carrying the sample Name (== the DIA-NN `Run` stem, the prolfqua join key),
  # the shared quant directory pointer, and that sample's [Factor] columns.
  # Registered by the SushiApp base class as a single multi-row dataset whose
  # parent is the DIANN child (so it is a grandchild of the raw dataset).
  def grandchild_datasets
    quant_dir = @params['enable_step_c'].to_s == 'true' ? 'out-DIANN_quantC' : 'out-DIANN_quantB'
    quant_path = File.join(@result_dir, quant_dir)
    rows = @dataset.is_a?(Array) ? @dataset : [@dataset].compact
    # Name the grandchild dataset after the DIANN result, not the first sample.
    @params['grandchildName'] = "#{@params['name']}_perSample"
    rows.map do |row|
      stripped = Hash[*row.map { |k, v| [k.to_s.gsub(/\[.+\]/, '').strip, v] }.flatten]
      grandchild = {
        'Name'               => stripped['Name'],
        'DIANN Quant [File]' => quant_path
      }
      # per-sample [Factor] columns (keep the tagged key), plus lineage columns
      @inherit_tags.each do |tag|
        row.each { |k, v| grandchild[k] = v if k.to_s.include?("[#{tag}]") }
      end
      row.each { |k, v| grandchild[k] = v if @inherit_columns.include?(k.to_s.gsub(/\[.+\]/, '').strip) }
      grandchild
    end
  end

  def commands
    # Python entry point: from ezpyz_diann.app import EzAppDiann
    # (run_PyApp resolves Ruby app name 'DIANN' to module ezpyz_diann + class EzAppDiann)
    run_PyApp('DIANN', pixi_enabled: true)  # Name must match [name] in 'ezpyz_[name]' format
  end
end

if __FILE__ == $0

end
