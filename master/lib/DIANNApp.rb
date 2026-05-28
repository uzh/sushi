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
    # Dataset must carry one of: 'RAW' (Hubert's earlier ezRun convention)
    # or 'Thermo RAW' (B-Fabric standard column). ezMethodDIANN picks
    # whichever is present, so we only require 'Name' here.
    @required_columns = ['Name', 'RAW']
    @required_params = ['name', 'paramsTemplate']

    # Compute resources
    @params['cores'] = '64'
    @params['cores', 'context']  = 'slurm'
    @params['ram']   = '32'
    @params['ram',   'context']  = 'slurm'
    @params['scratch'] = '20'
    @params['scratch', 'context'] = 'slurm'
    @params['node'] = ['fgcz-c-050']
    @params['node', 'context'] = 'slurm'

    # Params.yml selection: bundled template OR user-supplied override path.
    # The R method ezMethodDIANN reads these and writes work/params.yml.
    @params['paramsTemplate'] = ['default-DIA', 'default-DDA']
    @params['customParamsYml'] = ''
    @params['order_fasta'] = ''

    # Run identity
    @params['name'] = 'DIANN_v23'
    @params['mail'] = ''
    @modules = ['Dev/R']  # snakemake env is activated via @conda_env below
    @inherit_columns = []
  end

  def next_dataset
    {
      'Name' => @params['name'],
      'Protein Abundances [Link]' => File.join(@result_dir, 'qc_result/proteinAbundances.html'),
      'Sample Sizes [Link]'       => File.join(@result_dir, 'qc_result/QC_sampleSizeEstimation.html'),
      'qc_result [File]'          => File.join(@result_dir, 'qc_result'),
      'DIANN Quant [File]'        => File.join(@result_dir, 'DIANN_quantC')
    }.merge(extract_columns(colnames: @inherit_columns))
  end

  def commands
    run_RApp('EzAppDIANN', conda_env: 'gi_snakemake8.20.5')
  end
end

if __FILE__ == $0

end
