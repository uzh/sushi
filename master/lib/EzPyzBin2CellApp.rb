#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class EzPyzBin2CellApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'EzPyzBin2CellApp'
    @params['process_mode'] = 'SAMPLE'
    @analysis_category = 'Spatial'
    @description =<<-EOS
    A Bin2Cell app for VisiumHD data. Bin2cell is a tool for converting binned spatial transcriptomics data into single-cell format.<br/>
    The code is based on the original Bin2Cell app, find their github page here https://github.com/Teichlab/bin2cell?tab=readme-ov-file <br/>
    The Paper by Polanski et al is available here: https://academic.oup.com/bioinformatics/article/40/9/btae546/7754061?login=true <br/>
    This app is running fully on Python, therefore, check the implementation in ezPyz: https://github.com/fgcz/EzPyzApps .there you can find the required information and the implementation of the app.<br/>
    The compiled code lives in the conda environment `gi_bin2cell`<br/>
    EOS
    @required_columns = ['Name','BinnedOutputs2um','SourceImage','SpaceRanger']
    @required_params = ['name']
    @params['mpp'] = '0.5'
    @params['mpp','description'] = 'Microns per pixel. This is the pixel size of the image in microns.'
    @params['prob_thresh_HE'] = '0.01'
    @params['prob_thresh_HE','description'] = 'lowering `prob_thresh` to make the calls less stringent is recommended, for H&E image segmentation'
    @params['prob_thresh_GEX'] = '0.05'
    @params['prob_thresh_GEX','description'] = 'lowering `prob_thresh` to make the calls less stringent is recommended, for gene expression image segmentation'
    @params['roi_x'] = '-1'
    @params['roi_x','description'] = 'X coordinate of the region of interest (ROI) in pixels. If set to -1, the center of the image will be used.'
    @params['roi_x', "context"] = "EzPyzBin2CellApp"
    @params['roi_y'] = '-1'
    @params['roi_y','description'] = 'Y coordinate of the region of interest (ROI) in pixels. If set to -1, the center of the image will be used.'
    @params['roi_y', "context"] = "EzPyzBin2CellApp"
    @params['max_bin_distance'] = '2'
    @params['max_bin_distance','description'] = 'Distance (in bins) to expand the nuclei by'
    @params['max_bin_distance', "context"] = "EzPyzBin2CellApp"
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '30'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '100'
    @params['scratch', "context"] = "slurm"
    @params['name'] = 'Bin2CellApp'
    @params['mail'] = ""
    @modules = []
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def next_dataset    
    dir_name = "#{@params['name']}_#{@dataset['Name']}"
    report_dir = File.join(@result_dir, dir_name)
    {'Name'=> dir_name,
    'Bin2Cell [File]'=>report_dir,
    'Figures [Link]'=>File.join(report_dir, 'figures'),
    'Anndata [Link]'=>File.join(report_dir, 'cdata.h5ad'),
    'Stardist [Link]'=>File.join(report_dir, 'stardist')
    }
  end
  def commands
    run_PyApp("Bin2cell",conda_env: 'gi_bin2cell')
  end
end
