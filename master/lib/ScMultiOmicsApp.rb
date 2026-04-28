#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class ScMultiOmicsApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'ScMultiOmics'
    @params['process_mode'] = 'SAMPLE'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
Downstream multi-omics extension on top of an annotated ScSeurat scData.qs2.<br/>
Auto-discovers ADT, VDJ, and ATAC siblings of the original CountMatrix and
produces a multi-assay Seurat object plus an HTML report. RNA assay stays
default so existing exploreSC apps render unchanged.<br/>
Phase 1 supports ADT only; VDJ + ATAC + WNN land in subsequent phases.
    EOS
    # Required columns:
    # - CountMatrix: original CellRanger Multi/ARC count matrix (mtx dir or H5 parent).
    # - Report or SC Cluster Report: ScSeurat output dir containing scData.qs2.
    @required_columns = ['Name', 'Species', 'refBuild', 'CountMatrix']
    @required_params = ['name']
    # SLURM resources
    @params['cores'] = '4'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '64'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '50'
    @params['scratch', "context"] = "slurm"
    @params['name'] = 'ScMultiOmics'
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "reference genome assembly"
    # --- ADT options ---
    @params['adtNorm', 'hr-header'] = "ADT"
    @params['adtNorm'] = ['ADTnorm', 'CLR']
    @params['adtNorm', 'description'] = "ADT normalization method. ADTnorm is recommended for marker-rich panels; CLR is faster."
    @params['npcsADT'] = 18
    @params['npcsADT', 'description'] = "Number of PCs to compute on the ADT assay."
    # --- VDJ options (Phase 2) ---
    @params['vdjChain', 'hr-header'] = "VDJ (Phase 2)"
    @params['vdjChain'] = ['auto', 'TCR', 'BCR', 'both']
    @params['vdjChain', 'description'] = "Which VDJ chains to load. 'auto' detects siblings (vdj_t / vdj_b)."
    # --- WNN options (Phase 3) ---
    @params['runWNN', 'hr-header'] = "WNN (Phase 3)"
    @params['runWNN'] = true
    @params['runWNN', 'description'] = "Run WNN integration when 2+ dimensional modalities are present."
    # --- Misc ---
    @params['specialOptions'] = ''
    @params['mail'] = ''
    @modules = ["Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def preprocess
    @random_string = (1..12).map{[*('a'..'z')].sample}.join
  end
  def next_dataset
    report_file = File.join(@result_dir, "#{@dataset['Name']}_ScMultiOmicsReport")
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@dataset['Name'],
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'Static Report [Link]'=>report_link,
     'Report [File]'=>report_file,
     'ScMultiOmics [Link]'=>File.join(report_file, "scMultiData.qs2"),
    }.merge(extract_columns(@inherit_tags))
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
  end
  def commands
    run_RApp("EzAppScMultiOmics")
  end
end

if __FILE__ == $0

end
