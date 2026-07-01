#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class HUMAnNApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'HUMAnN'
    @params['process_mode'] = 'SAMPLE'
    @analysis_category = 'Metagenomics'
    @description =<<-EOS
HUMAnN 4 alpha functional profiling of shotgun metagenomes. Accepts either a
MetaPhlAn taxonomic profile (any DB) or a Bracken/Kraken report as upstream
input. HUMAnN mode is auto-detected from the upstream profile vintage:
<ul>
<li><b>Full HUMAnN</b> (tier-1 ChocoPhlAn + tier-2 UniRef90-EC; SGB-stratified output)
when the upstream MetaPhlAn profile carries header
<code>#mpa_vOct22_CHOCOPhlAnSGB_202403</code>.</li>
<li><b>Translated-only HUMAnN</b> (tier-2 only, community-level output) for any other
MetaPhlAn vintage, or for Bracken/Kraken input (the report is converted to MetaPhlAn-format
internally).</li>
</ul>
Outputs CPM-normalized gene families, MetaCyc reactions, MetaCyc pathways, an optional
KEGG KO regroup, plus a self-contained per-sample HTML report.
See <a href='https://github.com/biobakery/humann'>HUMAnN</a> and
<a href='https://forum.biobakery.org/t/metaphlan-4-humann-4-compatibility/8523'>HUMAnN 4 compatibility notes</a>.
EOS

    # XOR-of-AND: accept any of MetaPhlAnProfile, BrackenReport, KrakenReport.
    # Read2 is added by preprocess when paired=true.
    @required_columns = [
      ['Name', 'Read1', 'MetaPhlAnProfile'],
      ['Name', 'Read1', 'BrackenReport'],
      ['Name', 'Read1', 'KrakenReport']
    ]
    @required_params = ['paired']

    # ---- slurm --------------------------------------------------------------
    @params['cores']   = '32'
    @params['cores',   "context"] = "slurm"
    @params['ram']     = '96'
    @params['ram',     "context"] = "slurm"
    @params['scratch'] = '500'
    @params['scratch', "context"] = "slurm"

    # ---- HUMAnN core --------------------------------------------------------
    @params['paired'] = false
    @params['paired', "context"] = "HUMAnN"

    @params['forceTranslatedOnly'] = false
    @params['forceTranslatedOnly', 'description'] = 'Force translated-only mode. Default FALSE.'
    @params['forceTranslatedOnly', "context"] = "HUMAnN"

    @params['humannDbRoot'] = humann_db_choices
    @params['humannDbRoot', 'description'] = 'Directory containing chocophlan/, uniref/, utility_mapping/ subdirs for HUMAnN 4. Dropdown auto-populated from /srv/GT/databases/humann/ (only subdirs carrying all three required subdirs are listed; size shown after the name). CHOCOPhlAn vOct22 is currently the only database compatible with HUMAnN v4. The pangenomes in chocophlan/ subdirectory are used only when the input is a MetaPhlAn profile generated against CHOCOPhlAn vOct22; for any other upstream (newer MetaPhlAn database, Kraken/Bracken), HUMAnN runs in translated-only mode (using only uniref/ + utility_mapping/ subdirectories).'
    @params['humannDbRoot', "context"] = "HUMAnN"

    @params['keepStratifiedOutput'] = true
    @params['keepStratifiedOutput', 'description'] = 'Keep SGB-stratified rows in full-mode outputs (per-species functional contributions). Ignored in translated-only mode. Default TRUE.'
    @params['keepStratifiedOutput', "context"] = "HUMAnN"

    @params['regroupKEGGKO'] = true
    @params['regroupKEGGKO', 'description'] = 'Run UniRef90 -> KEGG KO regroup as an extra output table. HUMAnN 4 alpha utility_mapping only covers EC-annotated KOs, so coverage is partial. Default TRUE.'
    @params['regroupKEGGKO', "context"] = "HUMAnN"

    # ---- fastp trimming (copied from MetaPhlAnApp — same upstream cleanup) -
    @params['trimAdapter', 'hr'] = true
    @params['trimAdapter'] = true
    @params['trimAdapter', "context"] = "OpenGene/fastp"
    @params['trim_front1'] = '0'
    @params['trim_front1', "context"] = "OpenGene/fastp"
    @params['trim_tail1'] = '0'
    @params['trim_tail1', "context"] = "OpenGene/fastp"
    @params['cut_front'] = false
    @params['cut_front', "context"] = "OpenGene/fastp"
    @params['cut_front_window_size'] = '4'
    @params['cut_front_window_size', "context"] = "OpenGene/fastp"
    @params['cut_front_mean_quality'] = '20'
    @params['cut_front_mean_quality', "context"] = "OpenGene/fastp"
    @params['cut_tail'] = false
    @params['cut_tail', "context"] = "OpenGene/fastp"
    @params['cut_tail_window_size'] = '4'
    @params['cut_tail_window_size', "context"] = "OpenGene/fastp"
    @params['cut_tail_mean_quality'] = '20'
    @params['cut_tail_mean_quality', "context"] = "OpenGene/fastp"
    @params['cut_right'] = false
    @params['cut_right', "context"] = "OpenGene/fastp"
    @params['cut_right_window_size'] = '4'
    @params['cut_right_window_size', "context"] = "OpenGene/fastp"
    @params['cut_right_mean_quality'] = '20'
    @params['cut_right_mean_quality', "context"] = "OpenGene/fastp"
    @params['average_qual'] = '0'
    @params['average_qual', "context"] = "OpenGene/fastp"
    @params['max_len1'] = '0'
    @params['max_len1', "context"] = "OpenGene/fastp"
    @params['max_len2'] = '0'
    @params['max_len2', "context"] = "OpenGene/fastp"
    @params['poly_x_min_len'] = '10'
    @params['poly_x_min_len', "context"] = "OpenGene/fastp"
    @params['length_required'] = '18'
    @params['length_required', "context"] = "OpenGene/fastp"
    @params['cmdOptionsFastp'] = ''
    @params['cmdOptionsFastp', "context"] = "OpenGene/fastp"

    @params['mail'] = ""
    @modules = []   # conda-based; see commands block
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end

  def preprocess
    if @params['paired']
      @required_columns.each { |row| row << 'Read2' unless row.include?('Read2') }
    end
  end

  def next_dataset
    name = @dataset['Name']
    # The R worker writes "<name>.mode.txt" (containing "full" or
    # "translated") into the scratch CWD alongside the TSVs, so we can
    # publish the detected mode as a dataset column without re-parsing
    # the HUMAnN log. At next_dataset() time we're on the compute node
    # with CWD = scratch, so read by basename. Fall back to empty if
    # the file is absent (e.g. an older worker that didn't write it).
    mode_basename = "#{name}.mode.txt"
    mode_val = File.exist?(mode_basename) ? File.read(mode_basename).strip : ''

    {
      'Name'                    => name,
      'GeneFamiliesCPM [File]'  => File.join(@result_dir, "#{name}_genefamilies_cpm.tsv"),
      'KEGGKO_CPM [File]'       => File.join(@result_dir, "#{name}_genefamilies_ko_cpm.tsv"),
      'ReactionsCPM [File]'     => File.join(@result_dir, "#{name}_reactions_cpm.tsv"),
      'PathAbundance [File]'    => File.join(@result_dir, "#{name}_pathabundance.tsv"),
      'Mode [Characteristic]'   => mode_val,
      # [File] triggers a g-req copy of 00index.html to gstore; [Link]
      # surfaces it as a clickable URL in the SUSHI web UI. Both are
      # needed — the [Link] alone would leave the HTML uncopied.
      'Static Report [File]'    => File.join(@result_dir, '00index.html'),
      'Static Report [Link]'    => File.join(@result_dir, '00index.html'),
      'Live Report [Link]'      => "http://fgcz-shiny.uzh.ch/exploreMetaTax?data=#{@result_dir}",
    }.merge(extract_columns(@inherit_tags))
  end

  def commands
    command  = "set +e\n"
    command << "source /misc/ngseq12/miniforge3/etc/profile.d/conda.sh\n"
    command << "conda activate gi_qiime2-amplicon-2026.4\n"
    command << "set -e\n"
    command << run_RApp("EzAppHUMAnN")
  end

  # Scan /srv/GT/databases/humann/ for subdirectories that look like a
  # valid HUMAnN reference root (must contain chocophlan/, uniref/ and
  # utility_mapping/ subdirs). Returns Array of [label, value] pairs;
  # value is the full path; label includes the basename + on-disk size.
  # Returns [] when nothing valid is found (the multi-select widget will
  # then be empty and SUSHI will refuse to submit).
  def humann_db_choices
    root = '/srv/GT/databases/humann'
    return [] unless Dir.exist?(root)
    required_subs = %w[chocophlan uniref utility_mapping]

    candidates = Dir.entries(root)
                    .reject { |e| e.start_with?('.') }
                    .map    { |e| File.join(root, e) }
                    .select { |p| File.directory?(p) }
                    .select { |p| required_subs.all? { |s| Dir.exist?(File.join(p, s)) } }
                    .sort

    candidates.map do |path|
      name  = File.basename(path)
      bytes = `du -sb "#{path}" 2>/dev/null`.to_i
      gb    = bytes.to_f / (1024**3)
      label = bytes > 0 ? "#{name} (#{format('%.1f', gb)} GB)" : name
      [label, path]
    end
  end
end

if __FILE__ == $0
end
