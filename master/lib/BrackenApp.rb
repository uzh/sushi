#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class BrackenApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Bracken'
    @analysis_category = 'Metagenomics'
    @description =<<-EOS
Bracken (Bayesian Reestimation of Abundance with KrakEN): re-estimates species/genus-level
abundance from an existing Kraken2 report using read-length-specific kmer distributions.
<a href='https://github.com/jenniferlu717/Bracken'>https://github.com/jenniferlu717/Bracken</a>
EOS

    @required_columns = ['Name','KrakenReport']
    @required_params = ['brackenDBOpt']
    # slurm
    @params['cores'] = '2'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '8'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '50'
    @params['scratch', "context"] = "slurm"

    @params['brackenDBOpt'] = bracken_db_choices
    @params['brackenDBOpt', 'description'] = 'Bracken DB + read length (one entry per databaseNmers.kmer_distrib file found under /srv/GT/databases/kraken2/). Format: <DBname>/<N>mers. The DB must match the one used to build the input Kraken report. The read length <N> should match the sequencing dataset read length (check FastQC/MultiQC).'
    @params['brackenDBOpt', "context"] = "Bracken"
    @params['brackenLevel'] = ['S', 'S1', 'G', 'F', 'O', 'C', 'P', 'D']
    @params['brackenLevel', 'description'] = 'Taxonomic level for abundance re-estimation (Bracken -l): D=Domain, P=Phylum, C=Class, O=Order, F=Family, G=Genus, S=Species, S1=strain.'
    @params['brackenLevel', "context"] = "Bracken"
    @params['brackenThreshold'] = '0'
    @params['brackenThreshold', 'description'] = 'Minimum number of reads required at the chosen level in the Kraken report before re-estimation (Bracken -t). Default 0 (no filtering).'
    @params['brackenThreshold', "context"] = "Bracken"
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'extra commandline options for bracken; do not specify any option already covered by the dedicated input fields'
    @params['cmdOptions', "context"] = "Bracken"
    @params['mail'] = ""
    @modules = ["Dev/R", "Tools/Bracken"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end

  def next_dataset
    name = @dataset['Name']
    # KrakenReport and Read1/Read2 are pass-through pointers to upstream
    # data — NOT produced by this Bracken run. Expose as [Link] (no copy)
    # so downstream apps like HUMAnN can still see them, but SUSHI's
    # job_footer (sushiApp.rb:611-630) doesn't try to g-req copy them
    # back into their own source folders (which would fail with
    # "destination path already exists" and kill the whole job because
    # the wrapper runs under set -e). Only files actually produced in
    # scratch by this run go under [File].
    out = {'Name'=>name,
     'BrackenAbundance [File]'=>File.join(@result_dir, "#{name}.bracken"),
     'BrackenReport [File]'=>File.join(@result_dir, "#{name}.bracken.report.txt"),
     'KrakenReport [Link]'=>@dataset['KrakenReport'],
     'Live Report [Link]'=>"http://fgcz-shiny.uzh.ch/exploreMetaTax?data=#{@result_dir}",
    }
    out['Read1 [Link]'] = @dataset['Read1'] if @dataset['Read1']
    out['Read2 [Link]'] = @dataset['Read2'] if @dataset['Read2']
    out.merge(extract_columns(@inherit_tags))
  end

  def commands
    run_RApp("EzAppBracken")
  end

  # One entry per databaseNmers.kmer_distrib file under
  # /srv/GT/databases/kraken2/, formatted as "<DBname>/<N>mers". Grouped by DB
  # (alphabetical, 'Standard' first when present); within a DB, sorted by N.
  # This guarantees every option corresponds to a real on-disk index file —
  # no invalid (DB, read-length) combination can be selected from the UI.
  def bracken_db_choices
    root = '/srv/GT/databases/kraken2'
    return [] unless Dir.exist?(root)

    pairs = Dir.glob(File.join(root, '*', 'database*mers.kmer_distrib'))
               .map do |f|
                 db = File.basename(File.dirname(f))
                 n  = File.basename(f)[/database(\d+)mers\.kmer_distrib/, 1]
                 [db, n.to_i]
               end
               .compact

    dbs = pairs.map(&:first).uniq.sort
    preferred = 'Standard'
    dbs = [preferred] + (dbs - [preferred]) if dbs.include?(preferred)

    dbs.flat_map do |db|
      pairs.select { |d, _| d == db }
           .map(&:last).sort
           .map { |n| "#{db}/#{n}mers" }
    end
  end
end

if __FILE__ == $0
end
