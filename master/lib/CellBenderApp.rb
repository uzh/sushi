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
    @params['cores', 'description'] = 'CPU threads. Only used when gpu=0 (passed to cellbender as --cpu-threads). With gpu=1, CellBender does inference on the GPU and this is mostly unused.'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '30'
    @params['ram', 'description'] = 'Host RAM in GB (not GPU VRAM). 30 GB is fine for single-modality scRNA-seq. For multiome / ARC inputs (GEX + ATAC peaks, 100k+ features) bump to 60-100 GB — otherwise the cgroup OOM-kills the job.'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '100'
    @params['scratch', 'description'] = 'Local /scratch in GB for CellBender checkpoints (ckpt.tar.gz, written every 2 epochs) and the output h5. 100 GB is plenty.'
    @params['scratch', "context"] = "slurm"
    @params['gpu'] = '0'
    @params['gpu', 'description'] = 'Set to 1 to run cellbender with --cuda on a GPU (much faster, ~10-40 min vs hours on CPU). Requires partition=GPU. If you set gpu=1, do NOT change partition.'
    @params['gpu', "context"] = "CellBender"
    @params['gpu_feature'] = 'L40S'
    @params['gpu_feature', 'description'] = 'Restrict to GPUs matching this SLURM feature tag (sbatch -C). Keep L40S: the conda env gi_cellbender_0.3.2 uses an old PyTorch that does NOT support the Blackwell GPU on fgcz-c-056 (sm_120). Leaving c-056 open would cause silent fallback to CPU and a numel overflow on large matrices. Leave empty only if you know the env is compatible with every GPU in the partition.'
    @params['gpu_feature', "context"] = "slurm"
    @params['name'] = 'CellBender'
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'Extra cellbender remove-background flags. Common tweaks: --expected-cells N (override cell prior), --total-droplets-included N (extend the empty-droplet tail), --fpr 0.01 (target FDR), --epochs 150. Leave empty to use defaults. Do NOT re-pass --input, --output, --cuda or --cpu-threads; they are set by the app.'
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
