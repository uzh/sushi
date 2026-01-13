#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class CellRangerMultiApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'CellRangerMulti'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
This wrapper runs <a href='https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/using/multi',>cellranger multi</a> in Single-library analysis mode.<br><br>

Note: that running this app usually requires manual curation of the input dataset. Column names for the corresponding 10X library are as follows.<br>
<table style="width: 35%;">
<tbody>
<tr><td>VDJ-T</td><td>VdjTDataDir [File]</td></tr>
<tr><td>VDJ-B</td><td>VdjBDataDir [File]</td></tr>
<tr><td>Feature Barcoding</td><td>FeatureDataDir [File]</td></tr>
<tr><td>Multiplexing</td><td>MultiDataDir [File]</td></tr>
</tbody>
</table>
<br> 
When specifying multiplexing, use our simple <a href='https://fgcz-shiny.uzh.ch/app/sample2barcode'>ShinyApp</a> can be used to provide the barcoding information with or without feature barcoding.                  
    EOS
    @required_columns = ['Name','RawDataDir','Species']
    @required_params = ['name', 'refBuild']
    @params['cores'] = ['8', '12', '16']
    @params['cores', "context"] = "slurm"
    @params['ram'] = ['60', '40', '80']
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '500'
    @params['scratch', "context"] = "slurm"
    @params['name'] = 'CellRangerMulti'
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "referfence genome assembly"
    @params['refFeatureFile'] = 'genes.gtf'
    @params['featureLevel'] = 'gene'
    @params['TenXLibrary'] = ['GEX', 'VDJ-T', 'VDJ-B', 'FeatureBarcoding', 'Multiplexing', 'fixedRNA']
    @params['TenXLibrary', 'description'] = "Which 10X libraries? Note: Not all library types can be processed simultaneously. See the <a href='https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/using/multi#when'>support page</a> for further details. E.g. for fixedRNA, must also specify GEX."
    @params['TenXLibrary', 'multi_selection'] = true
    @params['TenXLibrary', 'selected'] = ['GEX', 'Multiplexing']
    @params['MultiplexingType'] = {'select'=>'', 'On chip multiplexing (OCM)'=>'ocm', 'Hashing with Antibody Capture'=>'antibody'}
    @params['MultiplexingType', 'description'] = "Which type of 3' multiplexing technology is used? See the <a href='https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-3p-multi'>support page</a> for further details. For OCM (On-chip Multiplexing), no MultiplexBarcodeSet is needed - Sample2Barcode file must have 'ocm_barcode_ids' column (OB1-OB4). For Hashing, select the corresponding TotalSeq AntibodyCapture csv file."
    @params['MultiplexBarcodeSet'] = {'select'=>''}
    Dir["/srv/GT/databases/10x/CMO_files/*"].sort.select{|design| File.file?(design)}.each do |dir|
      @params['MultiplexBarcodeSet'][File.basename(dir)] = File.basename(dir)
    end
    @params['MultiplexBarcodeSet', 'description'] = 'Used when CellPlex libraries. New files needs to be installed under /srv/GT/databases/10x/CMO_files'
    @params['probesetFile'] =  {'select'=>''}
    Dir["/srv/GT/databases/10x_Probesets/Chromium/*"].sort.select{|design| File.file?(design)}.each do |dir|
      @params['probesetFile'][File.basename(dir)] = File.basename(dir)
    end
    @params['probesetFile', 'description'] = 'set it only for probe-based single cell fixed RNA profiling (FRP)'
    @params['customProbesFile'] = ''
    @params['customProbesFile', 'file_upload'] = true
    @params['customProbesFile', 'description'] = 'Custom probeset CSV-file according to 10X specifications (https://tinyurl.com/10xProbeSetCSVFormat). ONLY for probe-based single cell fixed RNA profiling (FRP). Note that all genes listed must have a corresponding entry in secondRef or controlSeqs. Custom probes must have the same length as the probes in the reference file.'
    @params['FeatureBarcodeFile'] = ''
    @params['FeatureBarcodeFile', 'file_upload'] = true
    @params['FeatureBarcodeFile', 'description'] = '(e.g. for CITEseq)'
    @params['includeIntrons'] = true
    @params['includeIntrons', 'description'] = 'set to false to reproduce the default behavior in cell ranger v6 and earlier (NOTE: Ignored for fixedRNA)'
    @params['expectedCells'] = ''
    @params['expectedCells', 'description'] = 'Expected number of recovered cells. Leave this free to let cellranger estimate the cell number.'
    @params['transcriptTypes'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA', 'long_noncoding', 'short_noncoding', 'pseudogene']
    @params['transcriptTypes', 'multi_selection'] = true
    @params['transcriptTypes', 'selected'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA']
    @params['secondRef'] = ''
    @params['secondRef', 'description'] = 'full path to fasta file with e.g. viralGenes'
    @params['keepBam'] = true
    @params['keepBam', 'description'] = 'Keep bam file produced by CellRanger? Usually it is not neccessary for downstream analyses'
    @params['cmdOptions'] = ''
    @params['cmdOptions', "context"] = "CellRangerMulti"
    @params['cmdOptions', 'description'] = 'specify the commandline options for CellRanger (e.g. --include-introns for single nuclei data); do not specify any option that is already covered by the dedicated input fields'
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Tools/seqtk", "Dev/R", "Dev/Python", "Tools/samtools"]
    @params['CellRangerVersion'] = ["Aligner/CellRanger/9.0.0", "Aligner/CellRanger/10.0.0", "Aligner/CellRanger/8.0.1", "Aligner/CellRanger/7.1.0"]
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
      'transcriptTypes'=>@params['transcriptTypes'],
      'SCDataOrigin'=>'10X',
      'ResultDir [File,Link]'=>report_dir,
      'Read Count'=>@dataset['Read Count']
    }.merge(extract_columns(@inherit_tags))
    dataset
  end
  def commands
    command = "module load  #{@params["CellRangerVersion"]}\n"
    command << run_RApp("EzAppCellRangerMulti")
  end
end

if __FILE__ == $0

end
