#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class KrakenApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Kraken'
    @analysis_category = 'Metagenomics'
    @description =<<-EOS
Kraken taxonomic sequence classification system
<a href='https://ccb.jhu.edu/software/kraken2/index.shtml'>https://ccb.jhu.edu/software/kraken2/index.shtml</a>
EOS

    @required_columns = ['Name','Read1']
    @required_params = ['paired', 'krakenDBOpt']
    # optional params
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '100'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '500'
    @params['scratch', "context"] = "slurm"
    @params['paired'] = false
    @params['paired', "context"] = "Kraken"
    
    @params['krakenDBOpt'] = kraken_db_choices
    @params['krakenDBOpt', 'description'] = 'kraken database options with size on disk (auto-detected from /srv/GT/databases/kraken2/; see https://benlangmead.github.io/aws-indexes/k2). With multiDB checked, select multiple — k2 v2.17 will classify against all in one pass. WARNING: multi-DB mode is NCBI ONLY — every selected DB must share the NCBI taxonomy tree (custom DBs like GTDB, SILVA, UHGG cannot be mixed).'
    @params['krakenDBOpt', "context"] = "Kraken"
    @params['krakenDBOpt', 'multi_selection'] = true
    @params['multiDB'] = false
    @params['multiDB', 'description'] = 'enable multi-database classification (k2 v2.17, --db DB1,DB2,...). MULTI-DB IS NCBI ONLY — all selected DBs must share the NCBI taxonomy. If unchecked, only the first selected DB is used.'
    @params['multiDB', "context"] = "Kraken"
    @params['krakenConfidenceOpt'] = '0.0'
    @params['krakenConfidenceOpt', 'description'] = 'Confidence score threshold, between 0 and 1'
    @params['krakenConfidenceOpt', "context"] = "Kraken"
    @params['krakenPhredOpt'] = '0'
    @params['krakenPhredOpt', 'description'] = 'Phred score threshold'
    @params['krakenPhredOpt', "context"] = "Kraken"
    @params['minimum_hit_groups'] = ['2', '3', '4']
    @params['minimum_hit_groups', 'description'] = 'minimum number of hit groups (overlapping k-mers sharing the same minimizer) needed to make a call (kraken2 default: 2)'
    @params['minimum_hit_groups', "context"] = "Kraken"
    @params['report_minimizer_data'] = ['no', 'yes']
    @params['report_minimizer_data', 'description'] = 'add per-clade distinct-minimizer columns to the report (enables advanced filtering in exploreMetaTax)'
    @params['report_minimizer_data', "context"] = "Kraken"
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'specify other commandline options for kraken; do not specify any option that is already covered by the dedicated input fields'
    @params['cmdOptions', "context"] = "Kraken"
    
        # trimming options
    # general
    @params['trimAdapter', 'hr'] = true
    @params['trimAdapter'] = true
    @params['trimAdapter', "context"] = "OpenGene/fastp"
    # Fastp
    ## trimming
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
    @modules = ["QC/fastp", "Dev/R", "Tools/kraken/2.17.1", "Tools/Krona"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  end
  def next_dataset
    # Note: per-read .txt.gz and <name>_unclassified.fasta.gz are written to
    # the result dir but intentionally not surfaced here — only the aggregate
    # report and Krona are advertised downstream.
    {'Name'=>@dataset['Name'],
     'KronaReport [Link]'=>File.join(@result_dir, "#{@dataset['Name']}.html"),
     'KrakenReport [File]'=>File.join(@result_dir, "#{@dataset['Name']}.report.txt"),
     'KronaOutDir [File]'=>File.join(@result_dir, "#{@dataset['Name']}.html.files"),
     'KronaOut [File]'=>File.join(@result_dir, "#{@dataset['Name']}.html"),
     'Live Report [Link]'=>"http://fgcz-shiny.uzh.ch/exploreMetaTax?data=#{@result_dir}",
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    run_RApp("EzAppKraken")
  end

  # Scan /srv/GT/databases/kraken2/ for valid Kraken2 DB folders. A folder is
  # a valid DB if it contains 'hash.k2d'. Returns an ordered Hash
  # { "<folder> (<size> GB)" => "<folder>" } so the UI shows on-disk size while
  # the value passed to k2 is just the folder name. Standard is promoted first
  # when present; the rest are alphabetical. Returns {} if the root doesn't
  # exist or contains no valid DBs.
  def kraken_db_choices
    root = '/srv/GT/databases/kraken2'
    return {} unless Dir.exist?(root)

    dbs = Dir.entries(root)
             .reject { |e| e.start_with?('.') }
             .select { |e| File.exist?(File.join(root, e, 'hash.k2d')) }
             .sort

    preferred = 'Standard'
    dbs = [preferred] + (dbs - [preferred]) if dbs.include?(preferred)

    out = {}
    dbs.each do |name|
      bytes = `du -sb #{File.join(root, name)} 2>/dev/null`.to_i
      gb = (bytes.to_f / (1024**3))
      label = bytes > 0 ? "#{name} (#{format('%.1f', gb)} GB)" : name
      out[label] = name
    end
    out
  end
end

if __FILE__ == $0
end
