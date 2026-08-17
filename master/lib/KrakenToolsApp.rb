#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class KrakenToolsApp < SushiFabric::SushiApp
  # tool => process_mode + the input columns it needs on top of 'Name'.
  # Drives preprocess (process_mode, required_columns) and next_dataset.
  TOOL_SPEC = {
    'extract_kraken_reads' => {'mode' => 'SAMPLE',  'columns' => ['KrakenAssignments', 'KrakenReport', 'Read1']},
    'kreport2krona'        => {'mode' => 'SAMPLE',  'columns' => ['KrakenReport']},
    'kreport2mpa'          => {'mode' => 'SAMPLE',  'columns' => ['KrakenReport']},
    'combine_kreports'     => {'mode' => 'DATASET', 'columns' => ['KrakenReport']},
    'combine_mpa'          => {'mode' => 'DATASET', 'columns' => ['MpaProfile']},
    'filter_bracken.out'   => {'mode' => 'SAMPLE',  'columns' => ['BrackenAbundance']},
    'alpha_diversity'      => {'mode' => 'SAMPLE',  'columns' => ['BrackenAbundance']},
    'beta_diversity'       => {'mode' => 'DATASET', 'columns' => ['BrackenAbundance']},
    'make_kreport'         => {'mode' => 'SAMPLE',  'columns' => ['KrakenAssignments']},
    'make_ktaxonomy'       => {'mode' => 'DATASET', 'columns' => []},
    'fix_unmapped'         => {'mode' => 'DATASET', 'columns' => []},
  }

  def initialize
    super
    @name = 'KrakenTools'
    @citation = [
      'Lu J, Rincon N, Wood D E, Breitwieser F P, Pockrandt C, Langmead B, Salzberg S L, Steinegger M. Metagenome analysis using the Kraken software suite. Nature Protocols (2022). https://doi.org/10.1038/s41596-022-00738-y'
    ]
    @analysis_category = 'Metagenomics'
    @description =<<-EOS
KrakenTools 1.2.1: downstream analysis of Kraken2 / Bracken results.
<a href='https://github.com/jenniferlu717/KrakenTools'>https://github.com/jenniferlu717/KrakenTools</a>
<br/><br/>
Pick one script with <b>tool</b>. The parameter form shows the arguments of <i>all</i>
tools at once (SUSHI renders one flat form); only the block matching your chosen
tool is read, the rest is ignored. Each tool declares its own required input
columns, so submitting against the wrong upstream dataset fails immediately
with a message naming the missing column.
<br/><br/>
<b>Which dataset to run this on:</b>
<ul>
<li><code>extract_kraken_reads</code>, <code>make_kreport</code> &mdash; a Kraken dataset produced with <code>save_read_assignments = yes</code> (needs the per-read assignments).</li>
<li><code>kreport2krona</code>, <code>kreport2mpa</code>, <code>combine_kreports</code> &mdash; any Kraken dataset.</li>
<li><code>filter_bracken.out</code>, <code>alpha_diversity</code>, <code>beta_diversity</code> &mdash; a Bracken dataset.</li>
<li><code>combine_mpa</code> &mdash; the output of a previous <code>kreport2mpa</code> run of this app.</li>
<li><code>make_ktaxonomy</code>, <code>fix_unmapped</code> &mdash; database-build helpers: they read no dataset column, only the file paths you give below. Run them on any dataset of the project.</li>
</ul>
EOS

    @required_columns = ['Name']
    @required_params = ['tool']
    @params['process_mode'] = 'SAMPLE'
    @params['tool'] = TOOL_SPEC.keys
    @params['tool', 'description'] = 'which KrakenTools script to run. This also decides the process mode (per sample, or one job over the whole dataset) and the required input columns.'
    @params['cores'] = '4'
    @params['cores', 'context'] = 'slurm'
    @params['ram'] = '30'
    @params['ram', 'context'] = 'slurm'
    @params['scratch'] = '200'
    @params['scratch', 'context'] = 'slurm'

    @params['paired'] = false
    @params['paired', 'hr-header'] = 'extract_kraken_reads'
    @params['paired', 'description'] = 'extract_kraken_reads only: input is paired-end, so extract from Read1+Read2 and write two files.'
    @params['paired', 'context'] = 'extract_kraken_reads'
    @params['extract_taxids'] = ''
    @params['extract_taxids', 'description'] = 'REQUIRED for extract_kraken_reads. One or more NCBI taxids, comma- or space-separated (e.g. "9606 1280"). All of them are extracted in a single pass into one output file per sample.'
    @params['extract_taxids', 'context'] = 'extract_kraken_reads'
    @params['extract_includeChildren'] = ['yes', 'no']
    @params['extract_includeChildren', 'description'] = 'add --include-children: take the whole clade below each taxid. Almost always what you want &mdash; with "no", only reads assigned to that exact node are taken, which for an internal node is a small minority of its clade.'
    @params['extract_includeChildren', 'context'] = 'extract_kraken_reads'
    @params['extract_includeParents'] = ['no', 'yes']
    @params['extract_includeParents', 'description'] = 'add --include-parents: also take reads assigned to ancestors of each taxid.'
    @params['extract_includeParents', 'context'] = 'extract_kraken_reads'
    @params['extract_exclude'] = ['no', 'yes']
    @params['extract_exclude', 'description'] = 'add --exclude: DEPLETE the listed taxa instead of extracting them (host/contaminant removal). Note this depletes by k-mer assignment and will take genuine host reads with it, so it belongs in exploratory or microbial work, not upstream of a host count matrix.'
    @params['extract_exclude', 'context'] = 'extract_kraken_reads'
    @params['extract_outputFormat'] = ['fastq', 'fasta']
    @params['extract_outputFormat', 'description'] = 'fastq adds --fastq-output and the result is emitted as Read1/Read2 so it chains into any read-consuming app. fasta is the KrakenTools default and is emitted as ExtractedSeq instead.'
    @params['extract_outputFormat', 'context'] = 'extract_kraken_reads'
    @params['extract_max'] = '100000000'
    @params['extract_max', 'description'] = '--max: hard ceiling on reads written. This is a REAL cap, not "unlimited" &mdash; 100000000 is the KrakenTools default. Raise it for very deep libraries.'
    @params['extract_max', 'context'] = 'extract_kraken_reads'

    @params['k2krona_intermediateRanks'] = ['no', 'yes']
    @params['k2krona_intermediateRanks', 'hr-header'] = 'kreport2krona'
    @params['k2krona_intermediateRanks', 'description'] = '--intermediate-ranks: keep non-traditional ranks (those without a standard rank code).'
    @params['k2krona_intermediateRanks', 'context'] = 'kreport2krona'
    @params['k2krona_renderHtml'] = ['yes', 'no']
    @params['k2krona_renderHtml', 'description'] = 'also render the Krona text file to a self-contained HTML chart with ktImportText (module Tools/Krona). The .txt alone is only an intermediate.'
    @params['k2krona_renderHtml', 'context'] = 'kreport2krona'

    @params['k2mpa_displayHeader'] = ['no', 'yes']
    @params['k2mpa_displayHeader', 'hr-header'] = 'kreport2mpa'
    @params['k2mpa_displayHeader', 'description'] = '--display-header: write a header line naming the sample.'
    @params['k2mpa_displayHeader', 'context'] = 'kreport2mpa'
    @params['k2mpa_values'] = ['read_count', 'percentages']
    @params['k2mpa_values', 'description'] = 'report clade read counts (--read_count, default) or clade percentages (--percentages).'
    @params['k2mpa_values', 'context'] = 'kreport2mpa'
    @params['k2mpa_intermediateRanks'] = ['no', 'yes']
    @params['k2mpa_intermediateRanks', 'description'] = '--intermediate-ranks: include non-traditional ranks.'
    @params['k2mpa_intermediateRanks', 'context'] = 'kreport2mpa'
    @params['k2mpa_spaces'] = ['remove', 'keep']
    @params['k2mpa_spaces', 'description'] = 'remove (default, --remove-spaces) replaces spaces in taxon names with underscores; keep uses --keep-spaces.'
    @params['k2mpa_spaces', 'context'] = 'kreport2mpa'

    @params['combineK_displayHeaders'] = ['yes', 'no']
    @params['combineK_displayHeaders', 'hr-header'] = 'combine_kreports'
    @params['combineK_displayHeaders', 'description'] = 'yes = --display-headers (default), no = --no-headers.'
    @params['combineK_displayHeaders', 'context'] = 'combine_kreports'
    @params['combineK_onlyCombined'] = ['no', 'yes']
    @params['combineK_onlyCombined', 'description'] = '--only-combined: report only the summed columns, omitting the per-sample ones.'
    @params['combineK_onlyCombined', 'context'] = 'combine_kreports'

    @params['filterB_include'] = ''
    @params['filterB_include', 'hr-header'] = 'filter_bracken.out'
    @params['filterB_include', 'description'] = '--include: taxids to KEEP (comma- or space-separated). At least one of filterB_include / filterB_exclude must be set. Both may be used together, but no taxid may appear in both lists.'
    @params['filterB_include', 'context'] = 'filter_bracken.out'
    @params['filterB_exclude'] = ''
    @params['filterB_exclude', 'description'] = '--exclude: taxids to DROP (comma- or space-separated), with the abundances of the rest renormalised.'
    @params['filterB_exclude', 'context'] = 'filter_bracken.out'

    @params['alpha_metric'] = ['Sh', 'BP', 'F', 'Si', 'ISi']
    @params['alpha_metric', 'hr-header'] = 'alpha_diversity'
    @params['alpha_metric', 'description'] = '-a: Sh = Shannon (default), BP = Berger-Parker, F = Fisher, Si = Simpson, ISi = inverse Simpson. The script prints one value per sample; it is captured to a file.'
    @params['alpha_metric', 'context'] = 'alpha_diversity'

    @params['beta_type'] = ['bracken', 'single', 'simple', 'kreport', 'kreport2', 'krona']
    @params['beta_type', 'hr-header'] = 'beta_diversity'
    @params['beta_type', 'description'] = '--type: input file format. Use bracken for a Bracken dataset.'
    @params['beta_type', 'context'] = 'beta_diversity'
    @params['beta_level'] = ['all', 'S', 'G', 'F', 'O']
    @params['beta_level', 'description'] = '--level: taxonomic level to compare at (only used by the kreport/kreport2 types).'
    @params['beta_level', 'context'] = 'beta_diversity'
    @params['beta_cols'] = '1,2'
    @params['beta_cols', 'description'] = '--cols: 1-based name,count column pair, used only by --type single/simple.'
    @params['beta_cols', 'context'] = 'beta_diversity'

    @params['mkreport_taxonomy'] = ''
    @params['mkreport_taxonomy', 'hr-header'] = 'make_kreport'
    @params['mkreport_taxonomy', 'description'] = 'REQUIRED for make_kreport. Absolute path to a taxonomy file built by make_ktaxonomy (-t).'
    @params['mkreport_taxonomy', 'context'] = 'make_kreport'
    @params['mkreport_useReadLen'] = ['no', 'yes']
    @params['mkreport_useReadLen'] = ['no', 'yes']
    @params['mkreport_useReadLen', 'description'] = '--use-read-len: weight by summed read length instead of read count.'
    @params['mkreport_useReadLen', 'context'] = 'make_kreport'

    @params['ktax_nodes'] = ''
    @params['ktax_nodes', 'hr-header'] = 'make_ktaxonomy'
    @params['ktax_nodes', 'description'] = 'REQUIRED for make_ktaxonomy. Absolute path to the NCBI taxonomy nodes.dmp (--nodes), e.g. inside a Kraken2 DB at <db>/taxonomy/nodes.dmp.'
    @params['ktax_nodes', 'context'] = 'make_ktaxonomy'
    @params['ktax_names'] = ''
    @params['ktax_names', 'description'] = 'REQUIRED for make_ktaxonomy. Absolute path to names.dmp (--names).'
    @params['ktax_names', 'context'] = 'make_ktaxonomy'
    @params['ktax_seqid2taxid'] = ''
    @params['ktax_seqid2taxid', 'description'] = 'REQUIRED for make_ktaxonomy. Absolute path to the Kraken2 seqid2taxid.map (--seqid2taxid).'
    @params['ktax_seqid2taxid', 'context'] = 'make_ktaxonomy'

    @params['fixun_input'] = ''
    @params['fixun_input', 'hr-header'] = 'fix_unmapped'
    @params['fixun_input', 'description'] = 'REQUIRED for fix_unmapped. Absolute path to the unmapped-accession list from a failed kraken2-build (-i).'
    @params['fixun_input', 'context'] = 'fix_unmapped'
    @params['fixun_accession2taxid'] = ''
    @params['fixun_accession2taxid', 'description'] = 'REQUIRED for fix_unmapped. One or more absolute paths to accession2taxid files (--accession2taxid), space-separated.'
    @params['fixun_accession2taxid', 'context'] = 'fix_unmapped'

    @params['cmdOptions'] = ''
    @params['cmdOptions', 'hr-header'] = 'general'
    @params['cmdOptions', 'description'] = 'extra command line options appended to the chosen script; do not repeat anything already covered by a field above.'
    @params['cmdOptions', 'context'] = 'general'
    @params['mail'] = ''
    @modules = ['Dev/R', 'Tools/KrakenTools/1.2.1', 'Tools/Krona']
    @inherit_tags = ['Factor', 'B-Fabric', 'Characteristic']
  end

  def tool
    @params['tool'].to_s
  end

  def spec
    TOOL_SPEC[tool] or raise "KrakenToolsApp: unknown tool #{tool.inspect}"
  end

  # Derives process_mode and required_columns from the chosen tool, then validates
  # that tool's own parameters. Runs before check_required_columns (sushiApp.rb:1199
  # vs :1281), so a wrong tool/dataset pairing is rejected at submission time.
  def preprocess
    @params['process_mode'] = spec['mode']
    @required_columns = ['Name'] + spec['columns']
    @required_columns << 'Read2' if tool == 'extract_kraken_reads' and @params['paired']
    validate_params
  end

  def validate_params
    case tool
    when 'extract_kraken_reads'
      ids = @params['extract_taxids'].to_s.strip
      raise 'extract_kraken_reads: extract_taxids is empty; give at least one NCBI taxid.' if ids.empty?
      raise "extract_kraken_reads: extract_taxids must be numeric taxids separated by comma or space (got #{ids.inspect})." unless ids =~ /\A[\d[:space:],;]+\z/
      raise "extract_kraken_reads: extract_max must be a positive integer (got #{@params['extract_max'].inspect})." unless @params['extract_max'].to_s.strip =~ /\A[1-9]\d*\z/
    when 'filter_bracken.out'
      inc = @params['filterB_include'].to_s.strip
      exc = @params['filterB_exclude'].to_s.strip
      raise 'filter_bracken.out: set filterB_include or filterB_exclude (both empty would be a no-op).' if inc.empty? and exc.empty?
      [['filterB_include', inc], ['filterB_exclude', exc]].each do |key, val|
        next if val.empty?
        raise "filter_bracken.out: #{key} must be numeric taxids separated by comma or space (got #{val.inspect})." unless val =~ /\A[\d[:space:],;]+\z/
      end
      both = inc.split(/[\s,;]+/) & exc.split(/[\s,;]+/)
      raise "filter_bracken.out: taxid(s) #{both.join(', ')} appear in both filterB_include and filterB_exclude; the script rejects that." unless both.empty?
    when 'make_kreport'
      raise 'make_kreport: mkreport_taxonomy is required (the output of a make_ktaxonomy run).' if @params['mkreport_taxonomy'].to_s.strip.empty?
    when 'make_ktaxonomy'
      ['ktax_nodes', 'ktax_names', 'ktax_seqid2taxid'].each do |key|
        raise "make_ktaxonomy: #{key} is required." if @params[key].to_s.strip.empty?
      end
    when 'fix_unmapped'
      ['fixun_input', 'fixun_accession2taxid'].each do |key|
        raise "fix_unmapped: #{key} is required." if @params[key].to_s.strip.empty?
      end
    end
  end

  # Name of the single output row for the DATASET-mode tools, which produce one
  # file for the whole dataset rather than one per sample.
  def combined_name
    case tool
    when 'make_ktaxonomy' then 'ktaxonomy'
    when 'fix_unmapped'   then 'fix_unmapped'
    else 'combined'
    end
  end

  def extract_extension
    @params['extract_outputFormat'].to_s == 'fasta' ? 'fasta' : 'fastq'
  end

  def output_row_for(name)
    out = {'Name' => name}
    case tool
    when 'extract_kraken_reads'
      ext = extract_extension
      col = ext == 'fastq' ? 'Read' : 'ExtractedSeq'
      if @params['paired']
        out["#{col}1 [File]"] = File.join(@result_dir, "#{name}_extracted_R1.#{ext}.gz")
        out["#{col}2 [File]"] = File.join(@result_dir, "#{name}_extracted_R2.#{ext}.gz")
      else
        out["#{col}1 [File]"] = File.join(@result_dir, "#{name}_extracted.#{ext}.gz")
      end
      out['ExtractStats [File]'] = File.join(@result_dir, "#{name}.extract_stats.tsv")
    when 'kreport2krona'
      out['KronaText [File]'] = File.join(@result_dir, "#{name}.krona.txt")
      if @params['k2krona_renderHtml'].to_s == 'yes'
        out['KronaOut [File]'] = File.join(@result_dir, "#{name}.html")
        out['KronaReport [Link]'] = File.join(@result_dir, "#{name}.html")
      end
    when 'kreport2mpa'
      out['MpaProfile [File]'] = File.join(@result_dir, "#{name}.mpa.txt")
    when 'combine_kreports'
      out['CombinedKReport [File]'] = File.join(@result_dir, "#{name}.kreports_combined.txt")
    when 'combine_mpa'
      out['CombinedMpa [File]'] = File.join(@result_dir, "#{name}.mpa_combined.txt")
    when 'filter_bracken.out'
      out['BrackenAbundance [File]'] = File.join(@result_dir, "#{name}.filtered.bracken")
    when 'alpha_diversity'
      out['AlphaDiversity [File]'] = File.join(@result_dir, "#{name}.alpha_#{@params['alpha_metric']}.txt")
    when 'beta_diversity'
      out['BetaDiversity [File]'] = File.join(@result_dir, "#{name}.beta_diversity.txt")
    when 'make_kreport'
      out['KrakenReport [File]'] = File.join(@result_dir, "#{name}.report.txt")
    when 'make_ktaxonomy'
      out['KTaxonomy [File]'] = File.join(@result_dir, "#{name}.taxonomy.txt")
    when 'fix_unmapped'
      out['FixedUnmapped [File]'] = File.join(@result_dir, "#{name}.mapped.txt")
      out['StillUnmapped [File]'] = File.join(@result_dir, "#{name}.still_unmapped.txt")
    end
    out
  end

  def next_dataset
    if spec['mode'] == 'SAMPLE'
      name = @dataset['Name']
      output_row_for(name).merge(extract_columns(@inherit_tags, sample_name: name))
    else
      output_row_for(combined_name)
    end
  end

  def commands
    run_RApp('EzAppKrakenTools')
  end
end

if __FILE__ == $0
end
