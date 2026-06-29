#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class QIIME2App < SushiFabric::SushiApp
  def initialize
    super
    @name = 'QIIME2App'
    @analysis_category = 'Metagenomics'
    @description =<<-EOS
    Data processing with QIIME2 (amplicon-2026.4) for short-read Illumina 16S data.
    Includes DADA2 denoising, taxonomic classification against any reference
    database dropped under /srv/GT/databases/QIIME2/, alpha/beta diversity,
    optional fastp pre-trimming, optional PICRUSt2 functional prediction, and
    an ANCOM-BC-driven differential abundance Rmd report.
    <a href='https://qiime2.org/'>QIIME2 main website with tutorials and manuals.</a>

    Breaking change vs. the previous version: the symmetric `trim_left` and
    `truncate_len` parameters are replaced by `trim_left_f`/`trim_left_r` and
    `truncate_len_f`/`truncate_len_r`. The bash ANCOM step is dropped (ANCOM-BC
    now runs in the Rmd report).
EOS
    @params['process_mode'] = 'DATASET'
    @required_columns = ['Name', 'Read1']
    @required_params = ['paired','trim_left_f','truncate_len_f','sampling_depth','max_rarefaction_depth','min_freq','min_samples','group']

    # ---- slurm --------------------------------------------------------------
    @params['cores'] = '16'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '40'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '100'
    @params['scratch', "context"] = "slurm"

    # ---- paired / grouping --------------------------------------------------
    @params['paired'] = false
    @params['paired', 'description'] = 'whether the reads are paired end; if false then only Read1 is considered even if Read2 is available.'
    @params['paired', "context"] = "QIIME2"
    @params['group'] = true
    @params['group', 'description'] = 'There needs to be a group assignment column. Ensure the column name is in the format "NAME [Factor]"'
    @params['grouping'] = ''
    @params['grouping', 'description'] = 'Type in the name of the group assignment column. If the name is in the format "NAME [Factor]" type in the NAME (no [Factor] suffix).'

    # ---- DADA2 trimming (split: forward / reverse) --------------------------
    # Breaking change: legacy `trim_left` / `truncate_len` replaced by the
    # F/R pair below so paired-end runs can use asymmetric trimming.
    @params['trim_left_f'] = '0'
    @params['trim_left_f', 'description'] = "DADA2 trim-left for forward reads (5' bases to remove; usually = forward primer length when primers are not pre-trimmed)."
    @params['trim_left_r'] = '0'
    @params['trim_left_r', 'description'] = "DADA2 trim-left for reverse reads (5' bases to remove; usually = reverse primer length when primers are not pre-trimmed). Used only when paired=true."
    @params['truncate_len_f'] = '150'
    @params['truncate_len_f', 'description'] = 'DADA2 truncation length for forward reads. For single-end runs this is also the truncate length.'
    @params['truncate_len_r'] = '150'
    @params['truncate_len_r', 'description'] = 'DADA2 truncation length for reverse reads. Ignored when paired=false.'

    # ---- diversity / rarefaction --------------------------------------------
    @params['sampling_depth'] = '6000'
    @params['sampling_depth', 'description'] = 'Total frequency that each sample should be rarefied to prior to computing alpha/beta diversity metrics. Use feature-table summarize output to choose.'
    @params['max_rarefaction_depth'] = '4000'
    @params['max_rarefaction_depth', 'description'] = 'Maximum depth for generating alpha rarefaction curves.'

    # ---- differential abundance filtering -----------------------------------
    @params['min_freq'] = '1'
    @params['min_freq', 'description'] = 'Minimum total frequency a feature must reach to be retained for differential abundance calculation.'
    @params['min_samples'] = '1'
    @params['min_samples', 'description'] = 'Minimum number of samples a feature must appear in to be retained for differential abundance. Roughly 10% of samples is recommended.'

    # ---- reference database (auto-detected from /srv/GT/databases/QIIME2/) --
    # The dropdown is populated from on-disk <name>-seqs.qza / <name>-tax.qza
    # pairs. To register a new DB: drop both .qza files under
    # /srv/GT/databases/QIIME2/ — no code change needed.
    @params['database'] = qiime2_db_choices
    @params['database', 'description'] = 'Reference database (auto-detected). To add one, drop matching <name>-seqs.qza and <name>-tax.qza under /srv/GT/databases/QIIME2/.'

    # ---- primer / amplicon region -------------------------------------------
    # Primer region selector. forward_primer / reverse_primer arrays are kept
    # in the same order so Ctrl/Cmd-clicking the matching index in each
    # selects the conventional pair. V3-V4 promoted to first = default.
    @params['primer'] = ["V3-V4", "V1-V3(1)", "V1-V3(2)", "V3", "V4", "V4-V5", "V3-V5", "V4-V6"]
    @params['primer', 'description'] = '16S rRNA region used for training the Naive-Bayes classifier. Taxonomic classification improves when the classifier is trained on only the sequenced region.'
    @params['forward_primer'] = ["CCTACGGGNGGCWGCAG", "DAGAGTTTGATCMTGGCTCAG", "GAGAGTTTGATYMTGGCTCAG", "GATCCTACGGGAGGCAGCA", "GTGCCAGCMGCCGCGGTAA", "GTGCCAGCMGCCGCGGTAA", "CCTACGGGAGGCAGCAG", "GTGCCAGCMGCNGCGG3"]
    @params['forward_primer', 'description'] = 'Forward primer used for classifier extract-reads. Order tracks the `primer` dropdown above.'
    @params['reverse_primer'] = ["GACTACHVGGGTATCTAATCC", "TMTTACCGCGGCNGCTGGCAC", "ACCGCGGCTGCTGGCAC", "CTTACCGCGGCTGCTGGC", "GGACTACHVGGGTWTCTAAT", "CCGTCAATTCMTTTRAGTTT", "CCGTCAATTCMTTTRAGT", "GGGTTNCGNTCGTTG"]
    @params['reverse_primer', 'description'] = 'Reverse primer used for classifier extract-reads. Order tracks the `primer` dropdown above.'

    # ---- classifier extract-reads bounds + optional prebuilt ----------------
    @params['classifier_min_len'] = '350'
    @params['classifier_min_len', 'description'] = 'qiime feature-classifier extract-reads --p-min-length. Minimum length retained when carving the reference DB to the amplified region.'
    @params['classifier_max_len'] = '550'
    @params['classifier_max_len', 'description'] = 'qiime feature-classifier extract-reads --p-max-length. Upper bound on amplicon length kept for classifier training.'
    @params['classifier_path'] = ''
    @params['classifier_path', 'description'] = "Optional absolute path to a pre-built classifier .qza (e.g. under /srv/GT/databases/QIIME2/classifiers/). If set, skip extract-reads + fit-classifier-naive-bayes and use this artifact for classify-sklearn. Leave empty (or 'NONE') to train on the fly."

    # ---- optional pipeline steps (gate Rmd panels) --------------------------
    @params['run_fastp'] = true
    @params['run_fastp', 'description'] = 'Run fastp adapter trimming before DADA2 (paired-end uses --detect_adapter_for_pe). Quality and length filtering are disabled to leave DADA2 in charge.'
    @params['run_picrust2'] = false
    @params['run_picrust2', 'description'] = 'Run standalone picrust2_pipeline.py after taxonomy filtering. Produces ko_abundance.tsv, ec_abundance.tsv, pathway_abundance.tsv and gates the PICRUSt2 panel in the Rmd report.'

    # ---- output / notification ----------------------------------------------
    @params['name'] = 'QIIME2'
    @params['mail'] = ""
    @inherit_columns = ["Order Id"]
    # NOTE: no Lmod Dev/R loaded — the QIIME2 conda env (see commands) ships its
    # own R 4.5.2; pulling in system Dev/R prepended its bin to PATH and caused
    # data.table.so (built against conda R 4.5.2) to be loaded by system R 4.6
    # -> SETLENGTH undefined symbol. ezRun and all R deps must be installed
    # into the conda env's R library.
  end

  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
    # Backwards-compat: if a legacy YAML provides the old symmetric keys,
    # map them onto the new F params so old jobs deserialise cleanly.
    if @params['trim_left'] && !@params['trim_left'].to_s.empty?
      @params['trim_left_f'] = @params['trim_left']
      @params['trim_left_r'] = @params['trim_left']
    end
    if @params['truncate_len'] && !@params['truncate_len'].to_s.empty?
      @params['truncate_len_f'] = @params['truncate_len']
      @params['truncate_len_r'] = @params['truncate_len']
    end
  end

  def set_default_parameters
    @params['paired'] = dataset_has_column?('Read2')
  end

  def next_dataset
     nds = {'Name'=>@params['name']}
     nds['ResultDir [File]'] = File.join(@result_dir, 'Results_Folder/')
     nds['Static Report [Link]'] = File.join(@result_dir, 'Results_Folder/00index.html')
     nds['Demux Report [Link]'] = File.join(@result_dir, 'Results_Folder/demux_seqs.qzv.zip.folder/data/index.html')
     nds['Denoising stats [Link]'] = File.join(@result_dir, 'Results_Folder/dada2_denoising_stats.qzv.zip.folder/data/index.html')
     nds['Feature Table [Link]'] = File.join(@result_dir, 'Results_Folder/table-filtered.qzv.zip.folder/data/index.html')
     nds['Rep Seqs Report [Link]'] = File.join(@result_dir, 'Results_Folder/dada2_rep_set.qzv.zip.folder/data/index.html')
     nds['Taxonomy Barplot [Link]'] = File.join(@result_dir, 'Results_Folder/taxa-bar-plots.qzv.zip.folder/data/index.html')
     nds['Taxonomy List [Link]'] = File.join(@result_dir, 'Results_Folder/taxonomy.qzv.zip.folder/data/index.html')
     nds['Shannon Diversity [Link]'] = File.join(@result_dir, 'Results_Folder/shannon_group_significance.qzv.zip.folder/data/index.html')
     nds['Jaccard Diversity [Link]'] = File.join(@result_dir, 'Results_Folder/jaccard_group_significance.qzv.zip.folder/data/index.html')
     nds['Bray Curtis Diversity [Link]'] = File.join(@result_dir, 'Results_Folder/bray_curtis_group_significance.qzv.zip.folder/data/index.html')
     nds['Jaccard Emperor Plot [Link]'] = File.join(@result_dir, 'Results_Folder/jaccard_emperor_plot.qzv.zip.folder/data/index.html')
     nds['Bray Curtis Emperor Plot [Link]'] = File.join(@result_dir, 'Results_Folder/bray_curtis_emperor_plot.qzv.zip.folder/data/index.html')
     nds['Alpha rarefaction [Link]'] = File.join(@result_dir, 'Results_Folder/alpha-rarefaction.qzv.zip.folder/data/index.html')
     # Note: ANCOM bash step removed; differential abundance now runs in the Rmd
     # (ANCOM-BC). No ancom_group.qzv link is exposed here anymore.
     nds.merge(extract_columns(colnames: @inherit_columns))
  end

  def commands
     run_RApp("EzAppQIIME2", conda_env: "gi_qiime2-amplicon-2026.4")
  end

  # Scan /srv/GT/databases/QIIME2/ for <name>-seqs.qza files that have a
  # matching <name>-tax.qza alongside them. Returns the bare <name> labels
  # sorted alphabetically. Returns [] if the root doesn't exist or no valid
  # pair is present.
  def qiime2_db_choices
    root = '/srv/GT/databases/QIIME2'
    return [] unless Dir.exist?(root)
    Dir.entries(root)
       .select { |e| e.end_with?('-seqs.qza') }
       .map    { |e| e.sub(/-seqs\.qza\z/, '') }
       .select { |name| File.exist?(File.join(root, "#{name}-tax.qza")) }
       .sort
  end
end

if __FILE__ == $0

end
