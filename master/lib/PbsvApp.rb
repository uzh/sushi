#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class PbsvApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'pbsv'
    @analysis_category = 'Variants'
    @description =<<-EOS
pbsv - PacBio structural variant (SV) calling and analysis tools<br/>
<a href='https://github.com/PacificBiosciences/pbsv'>Github website</a>
EOS

    @required_columns = ['Name','BAM','BAI', 'refBuild']
    @required_params = ['refBuild','ReadOpt']
    # optional params
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '30'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '100'
    @params['scratch', "context"] = "slurm"
    @params['refBuild'] = ref_selector
    @params['refBuild', 'description'] = 'the genome refBuild and annotation to use as reference.'
    @params['refBuild', "context"] = "referfence genome assembly"
    @params['ReadOpt'] = 'HIFI'
    @params['ReadOpt', 'description'] = 'input read types: SUBREAD, CCS, HIFI. Default is HIFI'
    @params['ReadOpt', "context"] = "Pbsv"
    @params['region'] = ""
    @params['region', 'description'] = 'The region of the genome. You can give either a chromosome name or a region on a chromosome like chr1 or chr1:1000-2000'
    @params['region', "context"] = "Pbsv"
    @params['types'] = 'DEL,INS,INV,DUP,BND'
    @params['types', 'description'] = 'SV types to call: DEL,INS,INV,DUP,BND'
    @params['types', "context"] = "Pbsv"
    @params['minL'] = '20'
    @params['minL', 'description'] = 'Ignore SV with length shorter than the given length (bp)'
    @params['minL', "context"] = "Pbsv"
    @params['afOptions'] = '-q 20 -m 100'
    @params['afOptions', "context"] = "Pbsv"
    @params['afOptions', 'description'] = 'The options to filter alignments'
    @params['callOptions'] = '-A 3 -O 3 -S 1  -B 2 -P 20'
    @params['callOptions', 'description'] = 'The options to <a href=https://github.com/PacificBiosciences/pbsv> call</a>'
    @params['filterOptions'] = '--min-N-in-gap 50  --filter-near-reference-gap 1000 --filter-near-contig-end 1000'
    @params['filterOptions', 'description'] = 'The options to <a href=https://github.com/PacificBiosciences/pbsv>filter</a>'
    @params['cmdOptions'] = ""
    @params['cmdOptions', "context"] = "Pbsv"
    @params['mail'] = ""
    @modules = ["Variants/SURVIVOR", "Dev/R"]
  end
  def next_dataset
    {
     'Name'=>@dataset['Name'],
     'refBuild'=>@params['refBuild'],
     'OutDir [File]'=>File.join(@result_dir, "#{@dataset['Name']}"),
     'OutReport [Link]'=>File.join(@result_dir, "#{@dataset['Name']}", "#{@dataset['Name']}.html")
    }
  end
  def set_default_parameters
      @params['refBuild'] = @dataset[0]['refBuild']
  end
  def commands
    run_RApp("EzAppPbsv")
  end
end

if __FILE__ == $0

end

