#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class CellRangerApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'CellRangerCount'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
This wrapper runs <a href='https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count',>cellranger count</a> in Single-library analysis mode.
    EOS
    @required_columns = ['Name','RawDataDir','Species']
    @required_params = ['name', 'refBuild']
    @params['cores'] = '8'
    @params['ram'] = '40'
    @params['scratch'] = '200'
    @params['name'] = 'STARsolo'
    @params['refBuild'] = ref_selector
    @params['refFeatureFile'] = 'genes.gtf'
    @params['featureLevel'] = 'gene'
    #@params['TenXLibrary'] = ['GEX', 'VDJ']
    #@params['TenXLibrary', 'description'] = 'Which 10X library? GEX or VDJ.'
    #@params['scMode'] = ['SC', 'SN']
    #@params['scMode', 'description'] = 'Single-cell or single-nuclei?'
    @params['transcriptTypes'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA', 'long_noncoding', 'short_noncoding', 'pseudogene']
    @params['transcriptTypes', 'multi_selection'] = true
    @params['transcriptTypes', 'selected'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA']
    @params['controlSeqs'] = ''
    @params['controlSeqs', 'description'] = 'The extra control sequences (such as spikein sequences) available in https://fgcz-gstore.uzh.ch/reference/controlSeqs.fa'
    
    ## STARsolo parameters
    @params['soloType'] = ['CB_UMI_Simple','CB_UMI_Complex']
    @params['soloType', 'description'] = 'CB_UMI_Simple (a.k.a. Droplet), CB_UMI_Complex (e.g. Droplet).'
    @params['chemistry'] = ['SC3Pv3', 'SC3Pv2', 'SC3Pv1']
    @params['chemistry', 'description'] = 'Assay configuration and barcode whitelist.'
    @params['soloCBstart'] = '1'
    @params['soloCBstart', 'description'] = 'Cell barcode start base.'
    @params['soloCBlen'] = '16'
    @params['soloCBstart', 'description'] = 'Cell barcode length.'
    @params['soloUMIstart'] = '17'
    @params['soloUMIstart','description'] = 'UMI start base.'
    @params['soloUMIlen'] = 'auto'
    @params['soloUMIlen', 'description'] = 'UMI length. Select *auto* to use default UMI length based on chemistry. Specify the length otherwise.'
    
    #@params['soloFeatures'] = ['Gene','SJ','GeneFull','Transcript3p']
    #@params['soloFeatures', 'multi_selection'] = true
    #@params['soloFeatures','description'] = 'Genomic features for which the UMI counts per Cell Barcode are collected.'
    #@params['soloUMIfiltering'] = ['MultiGeneUMI','-']
    #@params['soloUMIfiltering','description'] = 'Type of UMI filtering. "MultiGeneUMI" removes lower-count UMIs that map to more than one gene; "-" performs basic filtering removing UMIs with N and homopolymers. '
    #@params['soloCBmatchWLtype'] = ['1MM_multi_pseudocounts','1MM_multi','1MM','Exact']
    #@params['soloCBmatchWLtype','description'] = 'matching the Cell Barcodes to the WhiteList.'
    
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'Specify the commandline options for CellRanger; do not specify any option that is already covered by the dedicated input fields'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    
    
    @modules = ["Dev/R"]
    @inherit_tags = ["B-Fabric"]
  end
  def set_default_parameters
  end
  def next_dataset
    report_dir = File.join(@result_dir,"#{@dataset['Name']}")
    if @params['TenXLibrary'] == "VDJ"
      dataset = {
        'Name'=>@dataset['Name'],
        'Species'=>@dataset['Species'],
        'refBuild'=>@params['refBuild'],
        'refFeatureFile'=>@params['refFeatureFile'],
        'featureLevel'=>@params['featureLevel'],
        'ResultDir [File]'=>report_dir,
        'Report [Link]'=>File.join(report_dir, 'web_summary.html')
      }.merge(extract_columns(@inherit_tags))
    else
      dataset = {
        'Name'=>@dataset['Name'],
        'Species'=>@dataset['Species'],
        'refBuild'=>@params['refBuild'],
        'refFeatureFile'=>@params['refFeatureFile'],
        'featureLevel'=>@params['featureLevel'],
        'transcriptTypes'=>@params['transcriptTypes'],
        'ResultDir [File]'=>report_dir,
        'Report [Link]'=>File.join(report_dir, 'web_summary.html'),
        'CountMatrix [Link]'=>File.join(report_dir, 'filtered_feature_bc_matrix')
      }.merge(extract_columns(@inherit_tags))
    end
    dataset
  end
  def commands
    run_RApp("EzAppCellRanger")
  end
end

if __FILE__ == $0

end
