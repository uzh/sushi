#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class MageckTestApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'MageckTest'
    # mageck test and MAGeCKFlute (which internally uses clusterProfiler)
    # unconditional. limma::alias2Symbol gated on species being hsa/mmu.
    # EnhancedVolcano deliberately excluded: no dedicated paper, pure plotting.
    @citation = [
      'Li, W. et al. MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. Genome Biology 15, 554 (2014). https://doi.org/10.1186/s13059-014-0554-4',
      'Wang, B. et al. Integrative analysis of pooled CRISPR genetic screens using MAGeCKFlute. Nature Protocols 14, 756-780 (2019). https://doi.org/10.1038/s41596-018-0113-7',
      'Wu, T. et al. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation 2(3), 100141 (2021). https://doi.org/10.1016/j.xinn.2021.100141',
      'Ritchie, M.E. et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47 (2015). https://doi.org/10.1093/nar/gkv007'
    ]
    @analysis_category = 'GenomeEditing'
    @description =<<-EOS
    Run test module in the tool Model-based Analysis of Genome-wide CRISPR-Cas9 Knockout (<a href='https://sourceforge.net/p/mageck/wiki/Home/'>MAGeCK</a>)
    Note: SampleNames starting with a number are not supported and will result in a crash of the software. 
    EOS
    @params['process_mode'] = 'DATASET'
    @required_columns = ['Name','Count', 'libName']
    @required_params = ['species','libName','sampleGroup','refGroup']
    # optional params
    @params['cores'] = ['1']
    @params['cores', "context"] = "slurm"
    @params['ram'] = ['8']
    @params['ram', "context"] = "slurm"
    @params['scratch'] = ['10']
    @params['scratch', "context"] = "slurm"
    @params['name'] = 'MAGeCK_Test'
    @params['species'] = ['hsa', 'mmu']
    @params['species', "context"] = "MageckTest"
    @params['libName'] = ''
    @params['libName', "context"] = "MageckTest"
    @params['specialOptions'] = ''
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'specify other commandline options for MAGeCK_Test; do not specify any option that is already covered by the dedicated input fields'
    @params['cmdOptions', "context"] = "MageckTest"
    @params['grouping'] = ''
    @params['sampleGroup'] = ''
    @params['sampleGroup', 'description'] = 'sampleGroup should be different from refGroup'
    @params['refGroup'] = ''
    @params['refGroup', 'description'] = 'refGroup should be different from sampleGroup'
    @params['useControls'] = true
    @params['useControls', 'description'] = 'use control sgRNAs for generating the null distribution of RRA' 
    @params['normalizationMethod'] = ['median', 'control', 'total', 'none']
    @params['geneLFCMethod'] = ['median','alphamedian','mean','alphamean','secondbest']
    @params['mail'] = ""
    @modules = ["Dev/R"]
    @inherit_columns = ["Order Id"]
  end
  def preprocess
  end
  def next_dataset
    @comparison = "#{@params['sampleGroup']}--over--#{@params['refGroup']}"
    @params['comparison'] = @comparison
    @params['name'] = @comparison
    report_folder = File.join(@result_dir, "#{@params['comparison']}")
    {'Name'=>@params['name'],
     'ReportData [File]'=>report_folder,
     'Report [Link]'=>File.join(report_folder, '00index.html')
    }.merge(extract_columns(colnames: @inherit_columns))
  end
  def set_default_parameters
    @params['libName'] = @dataset[0]['libName']
  end
  def commands
    run_RApp("EzAppMageckTest")
  end
end

if __FILE__ == $0

end

