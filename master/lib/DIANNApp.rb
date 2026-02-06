#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class DIANNApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'DIANN'
    @employee = true
    @analysis_category = 'Proteomics'
    @description =<<-EOS
Quality control after counting reads<br/>
    EOS
    @params['process_mode'] = 'DATASET'
    @required_columns = ['Name','RAW', 'Resource', 'Grouping Var']
    @required_params = []
    # optional params
    @params['cores'] = ['64']
    @params['cores', "context"] = "slurm"
    @params['ram'] = ['32']
    @params['ram', "context"] = "slurm"
    @params['scratch'] = ['20']
    @params['scratch', "context"] = "slurm"
    @params['name'] = 'DIANN_v23_DDA'
    @params['03_fasta_database_path'] = 'NONE'
    @params['03_fasta_use_custom'] = 'false'
    @params['order_fasta'] = "/srv/gstore/projects/p34486/o37485-order.fasta"
    @params['03b_additional_fasta_database_path'] = '/misc/fasta/p34486_Proteobench_TripleProteome_with_headers_20251008.fasta'
    @params['05_diann_is_dda'] = 'true'
    @params['06a_diann_mods_variable'] = '--var-mods 1 --var-mod UniMod:35,15.994915,M'
    @params['06b_diann_mods_no_peptidoforms'] = 'false'
    @params['06c_diann_mods_unimod4'] = 'true'
    @params['06d_diann_mods_met_excision'] = 'true'
    @params['07_diann_peptide_fragment_mz_max'] = '1800'
    @params['07_diann_peptide_fragment_mz_min'] = '200'
    @params['07_diann_peptide_max_length'] = '30'
    @params['07_diann_peptide_min_length'] = '6'
    @params['07_diann_peptide_precursor_charge_max'] = '3'
    @params['07_diann_peptide_precursor_charge_min'] = '2'
    @params['07_diann_peptide_precursor_mz_max'] = '1500'
    @params['07_diann_peptide_precursor_mz_min'] = '400'
    @params['08_diann_digestion_cut'] = 'K*,R*'
    @params['08_diann_digestion_missed_cleavages'] = '1'
    @params['09_diann_mass_acc_ms1'] = '15'
    @params['09_diann_mass_acc_ms2'] = '20'
    @params['10_diann_scoring_qvalue'] = '0.01'
    @params['11a_diann_protein_pg_level'] = 'protein_names_1'
    @params['11b_diann_protein_relaxed_prot_inf'] = 'true'
    @params['12a_diann_quantification_reanalyse'] = 'true'
    @params['12b_diann_quantification_no_norm'] = 'false'
    @params['13_diann_freestyle'] = 'None'
    @params['97_raw_converter'] = 'thermoraw'
    @params['98_diann_binary'] = 'diann-docker'
    @params['99_other_verbose'] = '1'
    @params['application_version'] = '2.3'    
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Dev/R", "Dev/uv"]
    @inherit_columns = []
  end
  def next_dataset
    report_file = File.join(@result_dir, @params['name'])
    #report_link = File.join(report_file, '00index.html')
    {'Name'=>@params['name'],
     'qc_result [File]'=> File.join(@result_dir, 'qc_result'),
     'Protein Abundances [Link]'=> File.join(@result_dir, 'qc_result/proteinAbundances.html'),
     'Sample Sizes [Link]'=> File.join(@result_dir, 'qc_result/QC_sampleSizeEstimation.html'),
     'Result [File]'=>File.join(@result_dir, "Result.zip"),
     'DIANN Quant [File]'=> File.join(@result_dir, "DIANN_quantB"),
    }.merge(extract_columns(colnames: @inherit_columns))
  end
  def commands
    run_RApp("EzAppDIANN", conda_env: "gi_snakemake8.20.5")
  end
end

if __FILE__ == $0

end

