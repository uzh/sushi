#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class DiffShotALDExApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'DiffShotALDEx'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Metagenomics'
    @description =<<-EOS
Differential abundance analysis on Bracken/Kraken2 shotgun metagenomics
profiles using **ALDEx2** (CLR + Monte-Carlo Dirichlet). With covariates the
GLM path (`aldex.glm`) is used; without covariates the classic t-test path.
ALDEx2 requires <b>integer read counts</b>, which is why this app only accepts
Bracken input — for MetaPhlAn-based DA use DiffShotMaAsLin3 or DiffShotLEfSe
(both work on the native MetaPhlAn relative-abundance output). Functional
(HUMAnN) input is not supported — use DiffShotMaAsLin3 for that.
Output is a single self-contained HTML report with volcano plots and result
tables at Species and Genus rank.
EOS

    @required_columns = ['Name', 'BrackenReport']
    @required_params  = ['grouping', 'sampleGroup', 'refGroup']

    # ---- slurm --------------------------------------------------------------
    @params['cores'] = ['4', '2', '8']
    @params['cores', "context"] = "slurm"
    @params['ram']   = ['32', '16', '64']
    @params['ram',   "context"] = "slurm"
    @params['scratch'] = ['50', '20', '100']
    @params['scratch', "context"] = "slurm"

    # ---- DA method (hardcoded, single-option dropdown) ---------------------
    @params['daMethod'] = ['ALDEx2']
    @params['daMethod', 'description'] = 'DA method run by this app. Hardcoded; use the DiffShotANCOMBC / DiffShotMaAsLin3 / DiffShotLEfSe apps for the other methods.'

    # ---- grouping (auto-populated from [Factor] columns by SUSHI) ----------
    @params['grouping'] = ''
    @params['grouping', 'description'] = 'Metadata [Factor] column to test on (the categorical variable of interest). Dropdown auto-populated from the input dataset.'
    @params['sampleGroup'] = ''
    @params['sampleGroup', 'description'] = 'Value of `grouping` treated as the contrast group. Must differ from refGroup.'
    @params['refGroup'] = ''
    @params['refGroup', 'description'] = 'Value of `grouping` used as the reference level. Coefficient signs read sampleGroup-vs-refGroup.'

    # ---- covariates (multi-select dropdown of remaining Factor cols) -------
    @params['daCovariates'] = ''
    @params['daCovariates', 'multi_selection'] = true
    @params['daCovariates', 'description'] = 'Multi-select dropdown auto-populated with every [Factor] column in the input dataset. Pick any subset to include as covariates; leaving empty runs the unadjusted t-test path.'

    # ---- filters -----------------------------------------------------------
    @params['prevalenceMin'] = 0.10
    @params['prevalenceMin', 'description'] = 'Prevalence filter: keep taxa seen (>1 read) in at least this fraction of samples. 0.10 = 10%.'
    @params['relabundMin'] = 0.10
    @params['relabundMin', 'description'] = 'Relative-abundance floor: keep taxa whose pooled read count is >= this percentage of the total. 0.10 = 0.10%.'

    # ---- thresholds shown on volcanos --------------------------------------
    @params['pValueThresh'] = 0.05
    @params['pValueThresh', 'description'] = 'p-value cut-off drawn as a horizontal line on volcano plots.'
    @params['log2RatioThresh'] = 1
    @params['log2RatioThresh', 'description'] = 'log2 fold-change cut-off drawn as vertical lines on volcano plots.'

    # ---- sample exclusion --------------------------------------------------
    @params['samplesToDrop'] = ''
    @params['samplesToDrop', 'description'] = 'Comma-separated list of SampleID values to drop entirely before any analysis (e.g. failed libraries).'

    # ---- kraken-biom -------------------------------------------------------
    @params['brackenMinRank'] = ['S', 'S1', 'G', 'F', 'O', 'C', 'P', 'D']
    @params['brackenMinRank', 'description'] = 'kraken-biom --min: lowest taxonomic rank collapsed to it. Should match the rank Bracken was run at (commonly S for species).'
    @params['brackenMaxRank'] = ['D', 'P', 'C', 'O', 'F', 'G', 'S']
    @params['brackenMaxRank', 'description'] = 'kraken-biom --max: highest taxonomic rank captured (D = Domain captures all reads).'

    # ---- misc --------------------------------------------------------------
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'Free-form extra options; reserved, currently unused.'
    @params['mail'] = ''
    @modules = []
    @inherit_tags = ['Factor', 'B-Fabric', 'Characteristic']
    @inherit_columns = ['Order Id']
  end

  def preprocess
    @random_string = (1..12).map{[*('a'..'z')].sample}.join
  end

  def next_dataset
    @comparison = "#{@params['sampleGroup']}--over--#{@params['refGroup']}"
    @params['comparison'] = @comparison
    @params['name']       = @comparison
    report_file = File.join(@result_dir, "#{@params['comparison']}")
    report_link = File.join(report_file, '00index.html')
    {'Name'                  => @comparison,
     'Static Report [Link]'  => report_link,
     'Report [File]'         => report_file,
    }.merge(extract_columns(colnames: @inherit_columns))
  end

  def set_default_parameters
    factor_cols = []
    if @dataset && @dataset[0]
      @dataset[0].each_key do |k|
        if k.to_s =~ /\[Factor\]\s*$/
          factor_cols << k.to_s.sub(/\s*\[Factor\]\s*$/, '').strip
        end
      end
    end
    @params['daCovariates'] = factor_cols
  end

  def commands
    command  = "set +e\n"
    command << "source /misc/ngseq12/miniforge3/etc/profile.d/conda.sh\n"
    command << "conda activate gi_qiime2-amplicon-2026.4\n"
    command << "set -e\n"
    command << run_RApp("EzAppDiffShot")
  end
end

if __FILE__ == $0
end
