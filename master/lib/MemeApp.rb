#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class MemeApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'MemeApp'
    @analysis_category = 'Motif'
    @description =<<-EOS
Perform motif discovery on DNA, RNA or protein datasets<br/>
<a href='http://meme-suite.org/tools/meme'>http://meme-suite.org/tools/meme</a>
    EOS
    @required_columns = ['Name','PeakSequences']
    @required_params = ['name']
    @params['cores'] = '1'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '10'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '20'
    @params['scratch', "context"] = "slurm"
    @params['name'] = 'MotifCheck_MEME'
    @params['motifDB'] = '-db /srv/GT/databases/MEME/motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_redundant_v2.meme -db /srv/GT/databases/MEME/motif_databases/MOUSE/uniprobe_mouse.meme -db /srv/GT/databases/MEME/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme'
    @params['motifDB', "context"] = "MemeApp"
    @params['cmdOptions'] = '-meme-mod zoops -minw 6 -maxw 30 -meme-nmotifs 5 -centrimo-score 5.0 -centrimo-ethresh 10.0'
    @params['cmdOptions', "context"] = "MemeApp"
    @params['filterPeaks'] = ['true', 'false']
    @params['filterPeaks', "context"] = "MemeApp"
    @params['specialOptions'] = 'distanceToTSS=2000 minFold=4'
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "referfence genome assembly"
    @params['mail'] = ''
    @modules = ['Dev/Perl','Dev/Python','Tools/MEME', "Dev/R", "Tools/BEDTools"]
  end
  def next_dataset
    meme_link = File.join(@result_dir, "#{@dataset['Name']}/#{@dataset['Name']}_meme-chip.html")
    {'Name'=>@dataset['Name'],
     'MEME Result [File]'=>File.join(@result_dir, "#{@dataset['Name']}"),
     'MEME Report [Link]'=>meme_link,
    }
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
  end
  def commands
    run_RApp("EzAppMEME")
  end
end

if __FILE__ == $0

end
