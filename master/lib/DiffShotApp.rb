#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class DiffShotApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'DiffShot'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Metagenomics'
    @description =<<-EOS
Differential abundance analysis on Bracken/Kraken2 shotgun metagenomics profiles.
Runs ALDEx2 (GLM), ANCOM-BC2, DESeq2, edgeR, MaAsLin3 and (optionally) LEfSe on
the same filtered BIOM table — every method except LEfSe can be covariate-adjusted.
Output is a single self-contained HTML report with volcano plots, overlap Venn,
and an interactive "Organisms of interest" search.
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

    # ---- grouping (auto-populated from [Factor] columns by SUSHI) ----------
    @params['grouping'] = ''
    @params['grouping', 'description'] = 'Metadata [Factor] column to test on (the categorical variable of interest). Dropdown auto-populated from the input dataset.'
    @params['sampleGroup'] = ''
    @params['sampleGroup', 'description'] = 'Value of `grouping` treated as the contrast group. Must differ from refGroup. Dropdown auto-populated from values in the chosen grouping column.'
    @params['refGroup'] = ''
    @params['refGroup', 'description'] = 'Value of `grouping` used as the reference level. Coefficient signs in every method read sampleGroup-vs-refGroup. Dropdown auto-populated from values in the chosen grouping column.'

    # ---- covariates (multi-select dropdown of remaining Factor cols) -------
    # Options are auto-populated from the input dataset in
    # set_default_parameters. The grouping column is filtered out at runtime
    # in R if the user accidentally selects it as a covariate too.
    @params['daCovariates'] = ''
    @params['daCovariates', 'multi_selection'] = true
    @params['daCovariates', 'description'] = 'Multi-select dropdown auto-populated with every [Factor] column in the input dataset. Pick any subset to include as covariates in ALDEx2 (GLM), ANCOM-BC2, DESeq2, edgeR and MaAsLin3. LEfSe ignores covariates by design. If you tick the same column you chose as grouping, it is silently dropped from the covariate list. Leave empty for the unadjusted model.'

    # ---- filters -----------------------------------------------------------
    @params['prevalenceMin'] = 0.10
    @params['prevalenceMin', 'description'] = 'Prevalence filter: keep taxa seen (>1 read) in at least this fraction of samples. 0.10 = 10%.'
    @params['relabundMin'] = 0.10
    @params['relabundMin', 'description'] = 'Relative-abundance floor: keep taxa whose pooled read count is >= this percentage of the total. 0.10 = 0.10%.'
    @params['libCut'] = 1000
    @params['libCut', 'description'] = 'ANCOM-BC2 minimum library size (lib_cut). Samples below this are excluded by ANCOM-BC2 only; other methods do not apply a sample-level floor here.'

    # ---- thresholds shown on volcanos --------------------------------------
    @params['pValueThresh'] = 0.05
    @params['pValueThresh', 'description'] = 'p-value cut-off drawn as a horizontal line on volcano plots and used in the overlap counts.'
    @params['log2RatioThresh'] = 1
    @params['log2RatioThresh', 'description'] = 'log2 fold-change cut-off drawn as vertical lines on volcano plots.'

    # ---- method-specific ---------------------------------------------------
    @params['maaslinMaxSig'] = 0.1
    @params['maaslinMaxSig', 'description'] = 'MaAsLin3 max_significance: q-value threshold above which associations are dropped from headline plots.'
    @params['runLefse'] = ['true', 'false']
    @params['runLefse', 'description'] = 'Whether to run LEfSe. Note: LEfSe does NOT accept covariates; when daCovariates is set the report still shows LEfSe but with a banner stating that its results are unadjusted.'

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
    # No Lmod 'Rversion' param: R itself comes from the activated conda env.
    @modules = []
    @inherit_tags = ['Factor', 'B-Fabric', 'Characteristic']
    @inherit_columns = ['Order Id']
  end

  def preprocess
    @comparison = "#{@params['sampleGroup']}--over--#{@params['refGroup']}"
    @params['comparison'] = @comparison
    @params['name']       = @comparison
  end

  def next_dataset
    report_file = File.join(@result_dir, @params['comparison'])
    report_link = File.join(report_file, '00index.html')
    {'Name'                  => @params['comparison'],
     'Static Report [Link]'  => report_link,
     'Report [File]'         => report_file,
    }.merge(extract_columns(colnames: @inherit_columns))
  end

  def set_default_parameters
    # Build the daCovariates multi-select options from every [Factor] column
    # in the input dataset header. Without this the multi-select renders empty.
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
    # R itself ships with the qiime2 conda env, so we activate that env BEFORE
    # invoking R. kraken-biom / biom-format / h5py are picked up the same way.
    command  = "source /misc/ngseq12/miniforge3/etc/profile.d/conda.sh\n"
    command << "conda activate gi_qiime2-amplicon-2026.4\n"
    command << run_RApp("EzAppDiffShot")
  end
end

if __FILE__ == $0
end
