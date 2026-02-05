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
When specifying multiplexing, use our simple <a href='https://fgcz-shiny.uzh.ch/app/sample2barcode'>ShinyApp</a> to provide the barcoding information with or without feature barcoding. See also the <a href='https://gitlab.bfabric.org/Genomics/paul-scripts/-/tree/main/.claude/skills/sample2barcode-generation'>Claude skill documentation</a> for detailed guidance.
    EOS
    @required_columns = ['Name','RawDataDir','Species']
    @required_params = ['name', 'refBuild']
    @params['cores'] = ['8', '12', '16']
    @params['cores', "context"] = "slurm"
    @params['ram'] = ['80', '40', '60', '100']
    @params['ram', "context"] = "slurm"
    @params['ram', 'description'] = "RAM per job in GB. Flex v2 with 96/384-plex multiplexing may require 100GB."
    @params['scratch'] = '500'
    @params['scratch', "context"] = "slurm"
    @params['name'] = 'CellRangerMulti'
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "referfence genome assembly"
    @params['refFeatureFile'] = 'genes.gtf'
    @params['featureLevel'] = 'gene'
    @params['TenXLibrary'] = ['GEX', 'VDJ-T', 'VDJ-B', 'FeatureBarcoding', 'Multiplexing', 'fixedRNA']
    @params['TenXLibrary', 'description'] = "Select 10X library types to process:<br>
• <b>GEX</b>: Gene Expression (required for most analyses)<br>
• <b>VDJ-T</b>: T-cell receptor (requires VdjTDataDir column)<br>
• <b>VDJ-B</b>: B-cell receptor (requires VdjBDataDir column)<br>
• <b>FeatureBarcoding</b>: CITE-seq/ADT (requires FeatureDataDir + FeatureBarcodeFile)<br>
• <b>Multiplexing</b>: HTO/CMO demux (HTO requires FeatureDataDir, CMO requires MultiDataDir, OCM needs neither)<br>
• <b>fixedRNA</b>: 10x Flex (Fixed RNA Profiling) - probe-based (requires GEX + probesetFile)<br>
<a href='https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-3p-multi'>10x Documentation</a>"
    @params['TenXLibrary', 'multi_selection'] = true
    @params['TenXLibrary', 'selected'] = ['GEX', 'Multiplexing']
    @params['MultiplexingType'] = {'select'=>'', 'OCM (On-Chip Multiplexing, max 4-plex)'=>'ocm', 'HTO (TotalSeq Hashtag Antibodies)'=>'antibody', 'CMO (10x CellPlex Lipid Barcodes)'=>'cmo'}
    @params['MultiplexingType', 'description'] = "Select 3' multiplexing technology:<br>
<b>OCM</b>: <code>ocm_barcode_ids</code> (OB1-OB4), no extra fastqs, no BarcodeSet<br>
<b>HTO</b>: <code>hashtag_ids</code> (e.g. B0251), needs FeatureDataDir + TotalSeq_AntibodyCapture<br>
<b>CMO</b>: <code>cmo_ids</code> (e.g. CMO301), needs MultiDataDir + 10x_CMO<br>
<a href='https://fgcz-shiny.uzh.ch/app/sample2barcode'>Generate Sample2Barcode</a>"
    @params['MultiplexBarcodeSet'] = {'select'=>''}
    Dir["/srv/GT/databases/10x/CMO_files/*"].sort.select{|design| File.file?(design)}.each do |dir|
      @params['MultiplexBarcodeSet'][File.basename(dir)] = File.basename(dir)
    end
    @params['MultiplexBarcodeSet', 'description'] = "Barcode reference file (match to MultiplexingType!):<br>
<b>HTO</b>: TotalSeqA/B/C_AntibodyCapture (A=universal, B=mouse, C=human)<br>
<b>CMO</b>: 10x_CMO (CMO301-CMO312)<br>
<b>OCM</b>: Leave empty"
    @params['probesetFile'] =  {'select'=>''}
    Dir["/srv/GT/databases/10x_Probesets/Chromium/*"].sort.select{|design| File.file?(design)}.each do |dir|
      @params['probesetFile'][File.basename(dir)] = File.basename(dir)
    end
    @params['probesetFile', 'description'] = "Required for 10x Flex (Fixed RNA Profiling). Available probe sets:<br>
• <b>v2.0.0</b>: Latest for Cell Ranger 10+ (Flex v2 chemistry)<br>
• <b>v1.1.x</b>: Updated coverage (MFRP/SFRP chemistry)<br>
• <b>v1.0.1</b>: Legacy probe set<br><br>
<b>For Flex multiplexing:</b> Sample2Barcode needs <code>probe_barcode_ids</code> column:<br>
• <b>Flex v1</b>: BC001-BC048 (e.g., BC001, BC002)<br>
• <b>Flex v2</b>: Plate-Well format (e.g., A-A01, A-B01) for 96/384-plex"
    @params['customProbesFile'] = ''
    @params['customProbesFile', 'file_upload'] = true
    @params['customProbesFile', 'description'] = 'Custom probeset CSV for 10x Flex (Fixed RNA). Format: <a href="https://tinyurl.com/10xProbeSetCSVFormat">10x specs</a>. Genes must have entries in secondRef or controlSeqs. Probes must match reference probe length.'
    @params['FeatureBarcodeFile'] = ''
    @params['FeatureBarcodeFile', 'file_upload'] = true
    @params['FeatureBarcodeFile', 'description'] = "Feature reference CSV for ADT protein quantification (CITE-seq). NOT for HTO demux (use MultiplexBarcodeSet instead).<br>
Columns: <code>id, name, read, pattern, sequence, feature_type</code> | <a href='https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-feature-ref-csv'>10x Format</a>"
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
    @params['CellRangerVersion'] = ["Aligner/CellRanger/10.0.0", "Aligner/CellRanger/9.0.0", "Aligner/CellRanger/8.0.1", "Aligner/CellRanger/7.1.0"]
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
