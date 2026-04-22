#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class CellBenderApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'CellBender'
    @params['process_mode'] = 'SAMPLE'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
CellBender is a software package for eliminating technical artifacts from high-throughput single-cell RNA sequencing (scRNA-seq) data. It is often also referred to by its French name, Le plieur de cellules. ToolLink: https://github.com/broadinstitute/CellBender<br/>
    EOS
    @required_columns = ['Name', 'Species', 'refBuild', 'refFeatureFile', 'CountMatrix']
    @required_params = ['name']
    # optional params
    @params['cores'] = '8'
    @params['cores', 'description'] = 'CPU threads. Not used for the CellBender inference itself (that runs on the GPU) but controls the SLURM cpu allocation that surrounds it.'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '30'
    @params['ram', 'description'] = 'Host RAM in GB (not GPU VRAM). 30 GB is fine for single-modality scRNA-seq. For multiome / ARC inputs (Peaks are dropped automatically before inference, so VRAM is fine) 30 GB host RAM remains sufficient.'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '100'
    @params['scratch', 'description'] = 'Local /scratch in GB for CellBender checkpoints (ckpt.tar.gz, written every 2 epochs) and the output h5. 100 GB is plenty.'
    @params['scratch', "context"] = "slurm"
    @params['partition'] = 'GPU_L40S'
    @params['partition', 'description'] = 'Pinned to GPU_L40S (fgcz-r-023 only). The conda env gi_cellbender_0.3.2 uses a PyTorch build that does NOT support the Blackwell GPU on fgcz-c-056 (sm_120), so the Blackwell partition would silently fall back to CPU and hang on large matrices. Do not change.'
    @params['partition', "context"] = "slurm"
    @params['gpu'] = '1'
    @params['gpu', 'description'] = 'Always 1 for this app. CellBender runs cellbender remove-background --cuda on the allocated L40S. Do not change.'
    @params['gpu', "context"] = "CellBender"
    @params['name'] = 'CellBender'
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'Extra cellbender remove-background flags. Common tweaks: --expected-cells N (override cell prior), --total-droplets-included N (extend the empty-droplet tail), --fpr 0.01 (target FDR), --epochs 150. Leave empty to use defaults. Do NOT re-pass --input, --output or --cuda; they are set by the app.'
    @params['cmdOptions', "context"] = "CellBender"
    @params['refBuild'] = ref_selector
    @params['refBuild', 'description'] = 'Reference genome assembly (informational only — CellBender itself does not align reads; the selection is stored as metadata in the output h5 and inherited by downstream apps).'
    @params['refBuild', "context"] = "referfence genome assembly"
    @params['refFeatureFile'] = 'genes.gtf'
    @params['refFeatureFile', 'description'] = 'Gene annotation file name, taken from the reference build (informational only — used as metadata, not for feature filtering).'
    @params['refFeatureFile', "context"] = "CellBender"
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def next_dataset
    report_file = File.join(@result_dir, "#{@dataset['Name']}")
    report_link = File.join(report_file, 'cellbender_report.html')
    {'Name'=>@dataset['Name'],
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'refFeatureFile'=>@params['refFeatureFile'],
     'Static Report [Link]'=>report_link,
     'ResultDir [File]'=>report_file,
     'CountMatrix [Link]'=>File.join(report_file,'cellbender_filtered_seurat.h5'),
     'UnfilteredCountMatrix [Link]'=>File.join(report_file,'cellbender_raw_seurat.h5'),
    }.merge(extract_columns(@inherit_tags))
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('refFeatureFile')
      @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    end
  end
  def commands
    run_RApp("EzAppCellBender",conda_env: "gi_cellbender_0.3.2")
  end
end

if __FILE__ == $0

end
