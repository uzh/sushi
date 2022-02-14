#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class SCOneSampleSCEApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'SCOneSampleSCE'
    @params['process_mode'] = 'SAMPLE'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
Single cell report<br/>
    EOS
    @required_columns = ['Name', 'Species', 'refBuild', 'CountMatrix', 'ResultDir']
    @required_params = ['name']
    # optional params
    @params['cores'] = '4'
    @params['ram'] = '20'
    @params['scratch'] = '50'
    @params['name'] = 'SCOneSampleSCE'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['vars.regress'] = ['none', 'CellCycle']
    @params['vars.regress', 'description'] = 'Variable to regress when modelling the gene variance. Choose CellCycle if it is a bias.'
    @params['tissue'] = []
    @params['tissue','multi_selection'] = true
    @params['tissue','all_selected'] = true
    @params['tissue', 'multi_selection_size'] = 10
    tissue = {}
    CSV.foreach("/srv/GT/databases/scGeneSets/all_cell_markers.txt", headers: true, col_sep: "\t") do |e|
      tissue[e["tissueType"]] = true
    end
    @params['tissue'] = tissue.keys.sort
    @params['tissue', 'description'] = 'Tissue the cells come from. Used in cell types identification for Human and Mouse organisms.'
    @params['resolution'] = [10, 20, 30, 40, 50]
    @params['resolution', 'description'] = 'A higher value will lead to more dense and less clusters.'
    @params['nreads'] = ''
    @params['nreads', 'description'] = 'Low quality cells have less than "nreads" reads. Only when applying fixed thresholds'
    @params['ngenes'] = ''
    @params['ngenes', 'description'] = 'Low quality cells have less than "ngenes" genes. Only when applying fixed thresholds'
    @params['perc_mito'] = ''
    @params['perc_mito', 'description'] = 'Low quality cells have more than "perc_mito" percent of mitochondrial genes. Only when applying fixed thresholds'
    @params['perc_ribo'] = ''
    @params['perc_ribo', 'description'] = 'Low quality cells have more than "perc_ribo" percent of ribosomal genes. Only when applying fixed thresholds'
    @params['cellsFraction'] = 0.05
    @params['cellsFraction', 'description'] = 'A gene will be kept if it is expressed in at least this fraction of cells'
    @params['nUMIs'] = 1
    @params['nUMIs', 'description'] = 'A gene will be kept if it has at least nUMIs in the fraction of cells specified before'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    #@params['Rversion'] = ["Dev/R/4.1.0", "Dev/R/4.0.4", "Dev/R/4.0.3"]
    @modules = ["Dev/R"]
  end
  def preprocess
    @random_string = (1..12).map{[*('a'..'z')].sample}.join
  end
  def next_dataset
    report_file = File.join(@result_dir, "#{@dataset['Name']}_SCReport")
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@dataset['Name'],
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'refFeatureFile'=>@params['refFeatureFile'],
     'Static Report [Link]'=>report_link,
     'SC Cluster Report [File]'=>report_file,
     'SC SCE [File]'=>File.join(report_file, "sce_h5"),
    }
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('refFeatureFile')
      @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    end
  end
  def commands
    #command = "module load #{@params["Rversion"]}\n"
    run_RApp("EzAppSCOneSampleSCE")
  end
end

if __FILE__ == $0

end

