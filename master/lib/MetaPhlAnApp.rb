#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class MetaPhlAnApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'MetaPhlAn'
    @analysis_category = 'Metagenomics'
    @description =<<-EOS
MetaPhlAn taxonomic profiling of microbial communities from metagenomic reads
<a href='https://github.com/biobakery/MetaPhlAn'>https://github.com/biobakery/MetaPhlAn</a>
EOS

    @required_columns = ['Name','Read1']
    @required_params = ['paired', 'metaphlanIndex']
    # optional params
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '64'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '200'
    @params['scratch', "context"] = "slurm"
    @params['paired'] = false
    @params['paired', "context"] = "MetaPhlAn"

    @params['metaphlanIndex'] = metaphlan_index_choices
    @params['metaphlanIndex', 'description'] = 'MetaPhlAn bowtie2 index basename (auto-detected from /srv/GT/databases/metaphlan_databases/, total on-disk size shown). The latest CHOCOPhlAnSGB DB pointed to by mpa_latest is promoted first when present. See <a href="http://segatalab.cibio.unitn.it/data/Database_links.html">MetaPhlAn DB index</a>.'
    @params['metaphlanIndex', "context"] = "MetaPhlAn"

    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'extra commandline options for metaphlan; do NOT specify --input_type, --db_dir, --index, --nproc, --mapout, -o, --tmp_dir (already set by the app).'
    @params['cmdOptions', "context"] = "MetaPhlAn"

    # trimming options (same as Kraken — fastp upstream)
    @params['trimAdapter', 'hr'] = true
    @params['trimAdapter'] = true
    @params['trimAdapter', "context"] = "OpenGene/fastp"
    @params['trim_front1'] = '0'
    @params['trim_front1','description'] = 'trimming how many bases in front for read1 (and read2), default is 0.'
    @params['trim_front1', "context"] = "OpenGene/fastp"
    @params['trim_tail1'] = '0'
    @params['trim_tail1','description'] = 'trimming how many bases in tail for read1 (and read2), default is 0.'
    @params['trim_tail1', "context"] = "OpenGene/fastp"
    @params['cut_front'] = false
    @params['cut_front','description'] = 'move a sliding window from front (5p) to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.'
    @params['cut_front', "context"] = "OpenGene/fastp"
    @params['cut_front_window_size'] = '4'
    @params['cut_front_window_size', "context"] = "OpenGene/fastp"
    @params['cut_front_mean_quality'] = '20'
    @params['cut_front_mean_quality', "context"] = "OpenGene/fastp"
    @params['cut_tail'] = false
    @params['cut_tail','description'] = 'move a sliding window from tail (3p) to front, drop the bases in the window if its mean quality < threshold, stop otherwise.'
    @params['cut_tail', "context"] = "OpenGene/fastp"
    @params['cut_tail_window_size'] = '4'
    @params['cut_tail_window_size', "context"] = "OpenGene/fastp"
    @params['cut_tail_mean_quality'] = '20'
    @params['cut_tail_mean_quality', "context"] = "OpenGene/fastp"
    @params['cut_right'] = false
    @params['cut_right','description'] = 'move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, and then stop.'
    @params['cut_right', "context"] = "OpenGene/fastp"
    @params['cut_right_window_size'] = '4'
    @params['cut_right_window_size', "context"] = "OpenGene/fastp"
    @params['cut_right_mean_quality'] = '20'
    @params['cut_right_mean_quality', "context"] = "OpenGene/fastp"
    @params['average_qual'] = '0'
    @params['average_qual','description'] = 'if one read average quality score <avg_qual>, then this read/pair is discarded. Default 0 means no requirement'
    @params['average_qual', "context"] = "OpenGene/fastp"
    @params['max_len1'] = '0'
    @params['max_len1','description'] = 'if read1 is longer than max_len1, then trim read1 at its tail to make it as long as max_len1. Default 0 means no limitation. If two reads are present, the same will apply to read2.'
    @params['max_len1', "context"] = "OpenGene/fastp"
    @params['max_len2'] = '0'
    @params['max_len2','description'] = 'if read1 is longer than max_len2, then trim read2 at its tail to make it as long as max_len1. Default 0 means no limitation.'
    @params['max_len2', "context"] = "OpenGene/fastp"
    @params['poly_x_min_len'] = '10'
    @params['poly_x_min_len','description'] = 'the minimum length to detect polyX in the read tail. 10 by default.'
    @params['poly_x_min_len', "context"] = "OpenGene/fastp"
    @params['length_required'] = '18'
    @params['length_required','description'] = 'reads shorter than length_required will be discarded.'
    @params['length_required', "context"] = "OpenGene/fastp"
    @params['cmdOptionsFastp'] = ''
    @params['cmdOptionsFastp', "context"] = "OpenGene/fastp"
    @params['mail'] = ""
    @modules = ["QC/fastp", "Dev/R", "Tools/MetaPhlAn/4.2.4"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end

  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  end

  def next_dataset
    # Output filename ends with _metaphlan.txt so exploreMetaTax auto-detects
    # it as a MetaPhlAn profile (see exploreMetaTax/app.R FORMAT_PATTERNS).
    {'Name'=>@dataset['Name'],
     'MetaPhlAnProfile [File]'=>File.join(@result_dir, "#{@dataset['Name']}_metaphlan.txt"),
     'BowtieMapout [File]'=>File.join(@result_dir, "#{@dataset['Name']}.bowtie2.bz2"),
     'Live Report [Link]'=>"http://fgcz-shiny.uzh.ch/exploreMetaTax?data=#{@result_dir}",
    }.merge(extract_columns(@inherit_tags))
  end

  def commands
    run_RApp("EzAppMetaPhlAn")
  end

  # Scan /srv/GT/databases/metaphlan_databases/ for bowtie2 index basenames.
  # Files lie flat in the folder; a basename is valid if <basename>.1.bt2(l)
  # AND <basename>.pkl both exist. Returns an Array of [label, value] pairs
  # (Array, not Hash — SUSHI's set_parameters.html.erb only renders the
  # multi-select widget for Arrays). The basename pointed to by mpa_latest
  # is promoted first when present; the rest are alphabetical. Returns []
  # if the root doesn't exist or contains no valid indexes.
  def metaphlan_index_choices
    root = '/srv/GT/databases/metaphlan_databases'
    return [] unless Dir.exist?(root)

    files = Dir.entries(root).reject { |e| e.start_with?('.') || File.directory?(File.join(root, e)) }
    basenames = files.flat_map { |f|
      m = f.match(/\A(.+?)\.(?:rev\.)?[1-4]\.bt2l?\z/)
      m ? [m[1]] : []
    }.uniq

    valid = basenames.select { |b|
      pkl = File.join(root, "#{b}.pkl")
      idx1 = ["#{b}.1.bt2l", "#{b}.1.bt2"].any? { |f| File.exist?(File.join(root, f)) }
      File.exist?(pkl) && idx1
    }.sort

    latest_file = File.join(root, 'mpa_latest')
    if File.exist?(latest_file)
      preferred = File.read(latest_file).strip
      valid = [preferred] + (valid - [preferred]) if valid.include?(preferred)
    end

    valid.map do |name|
      bytes = Dir.glob(File.join(root, "#{name}.*"))
                 .map { |p| File.file?(p) ? File.size(p) : 0 }.sum
      gb = (bytes.to_f / (1024**3))
      label = bytes > 0 ? "#{name} (#{format('%.1f', gb)} GB)" : name
      [label, name]
    end
  end
end

if __FILE__ == $0
end
