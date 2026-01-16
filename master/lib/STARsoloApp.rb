#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class STARsoloApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'STARsolo'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
This wrapper runs <a href='https://github.com/alexdobin/STAR/blob/2.7.3a/docs/STARsolo.md',>STARsolo</a> in Single-library analysis mode. Note that it only runs on Single Cell GEX 10X libraries.
    EOS
    ## general
    @required_columns = ['Name','RawDataDir','Species']
    @required_params = ['name', 'refBuild']
    @params['cores'] = '8'
    @params['cores', "context"] = "slurm"
    @params['ram'] = '60'
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '300'
    @params['scratch', "context"] = "slurm"
    @params['name'] = 'STARsolo'
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "referfence genome assembly"
    @params['refFeatureFile'] = 'genes.gtf'
    @params['refFeatureFile', "context"] = "STARsolo"
    @params['featureLevel'] = 'gene'
    @params['featureLevel', "context"] = "STARsolo"

    ## STARsolo parameters
    @params['soloType'] = ['CB_UMI_Simple','CB_UMI_Complex']
    @params['soloType', 'description'] = 'CB_UMI_Simple (a.k.a. Droplet), CB_UMI_Complex (e.g. Droplet).'
    @params['soloCBwhitelist'] = ['SC3Pv4', 'SC3Pv3', 'SC3Pv2', 'SC3Pv1']
    @params['soloCBwhitelist', 'description'] = 'Select barcode whitelist. For details on which one to use, see: https://kb.10xgenomics.com/s/article/115004506263-What-is-a-barcode-inclusion-list-formerly-barcode-whitelist'
    @params['soloCBstart'] = '1'
    @params['soloCBstart', 'description'] = 'Cell barcode start base.'
    @params['soloCBlen'] = '16'
    @params['soloCBstart', 'description'] = 'Cell barcode length.'
    @params['soloUMIstart'] = '17'
    @params['soloUMIstart','description'] = 'UMI start base.'
    @params['soloUMIlen'] = 'auto'
    @params['soloUMIlen', 'description'] = 'UMI length. Select *auto* to use default UMI length based on chemistry. Specify the length otherwise.'
    @params['soloFeatures'] = ['GeneFull_ExonOverIntron', 'GeneFull', 'GeneFull_Ex50pAS', 'SJ', 'Gene', 'Velocyto']
    @params['soloFeatures', 'multi_selection'] = true
    @params['soloFeatures', 'selected'] = ['GeneFull_ExonOverIntron', 'Gene', 'Velocyto']
    @params['soloFeatures', 'description'] = "Specify genomic features for which the UMI counts per Cell Barcode are collected. Note 'Velocyto' requires 'Gene' to be specified."
    @params['keepAlignment'] = false
    @params['transcriptTypes'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA', 'long_noncoding', 'short_noncoding', 'pseudogene']
    @params['transcriptTypes', 'multi_selection'] = true
    @params['transcriptTypes', 'selected'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA']
    @params['controlSeqs'] = ''
    @params['controlSeqs', 'description'] = 'The extra control sequences (such as spikein sequences) available in https://fgcz-gstore.uzh.ch/reference/controlSeqs.fa'
    
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'Specify the commandline options for CellRanger; do not specify any option that is already covered by the dedicated input fields'
    @params['cmdOptions', "context"] = "STARsolo"
    @params['specialOptions'] = ''
    @params['mail'] = ""
    
    
    @modules = ["Dev/R", "Tools/samtools", "Aligner/STAR","Aligner/CellRanger"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def set_default_parameters
  end
  def next_dataset
    report_dir = File.join(@result_dir,"#{@dataset['Name']}")
      dataset = {
        'Name'=>@dataset['Name'],
        'Species'=>@dataset['Species'],
        'refBuild'=>@params['refBuild'],
        'refFeatureFile'=>@params['refFeatureFile'],
        'featureLevel'=>@params['featureLevel'],
        'soloFeatures'=>@params['soloFeatures'],
        'transcriptTypes'=>@params['transcriptTypes'],
        'ResultDir [File]'=>report_dir,
        'CountMatrix [Link]'=>File.join(report_dir, 'filtered_feature_bc_matrix'),
        'UnfilteredCountMatrix [Link]'=>File.join(report_dir, 'raw_feature_bc_matrix')
      }.merge(extract_columns(@inherit_tags))
    dataset
  end
  def commands
    run_RApp("EzAppSTARsolo")
  end
end

if __FILE__ == $0

end
