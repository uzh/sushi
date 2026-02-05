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
    @params['ram'] = ['60', '80', '100', '40']
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
    @params['MultiplexingType', 'description'] = "Select the 3' multiplexing technology used:<br><br>
<b>OCM (On-Chip Multiplexing / Overhang)</b><br>
• Max 4 samples per GEM well (OB1-OB4)<br>
• Barcodes embedded in GEX reads - no separate fastqs needed<br>
• Sample2Barcode: <code>ocm_barcode_ids</code> column (unquoted, values: OB1-OB4)<br>
• Leave MultiplexBarcodeSet empty<br><br>
<b>HTO (TotalSeq Hashtag Antibodies)</b><br>
• Uses TotalSeq-A/B/C antibodies<br>
• Requires Antibody Capture fastqs: add <code>FeatureDataDir</code> column<br>
• Sample2Barcode: <code>hashtag_ids</code> column with TotalSeq IDs (e.g., \"B0251\")<br>
• Select TotalSeq_AntibodyCapture file in MultiplexBarcodeSet<br><br>
<b>CMO (10x CellPlex)</b><br>
• Uses 10x lipid-conjugated barcodes (CMO301-CMO312)<br>
• Requires Multiplexing Capture fastqs: add <code>MultiDataDir</code> column<br>
• Sample2Barcode: <code>cmo_ids</code> column with CMO IDs (e.g., \"CMO301\")<br>
• Select 10x_CMO file in MultiplexBarcodeSet<br><br>
Generate Sample2Barcode files: <a href='https://fgcz-shiny.uzh.ch/app/sample2barcode'>ShinyApp</a> | <a href='https://gitlab.bfabric.org/Genomics/paul-scripts/-/tree/main/.claude/skills/sample2barcode-generation'>Documentation</a>"
    @params['MultiplexBarcodeSet'] = {'select'=>''}
    Dir["/srv/GT/databases/10x/CMO_files/*"].sort.select{|design| File.file?(design)}.each do |dir|
      @params['MultiplexBarcodeSet'][File.basename(dir)] = File.basename(dir)
    end
    @params['MultiplexBarcodeSet', 'description'] = "Select the barcode reference file for your multiplexing technology:<br><br>
<b>⚠️ IMPORTANT: Match file to multiplexing method!</b><br><br>
<b>For HTO (Hashtag) demultiplexing - use *_AntibodyCapture.csv files:</b><br>
• <b>TotalSeqA_AntibodyCapture</b>: Universal hashtags<br>
• <b>TotalSeqB_AntibodyCapture</b>: Mouse-specific hashtags<br>
• <b>TotalSeqC_AntibodyCapture</b>: Human-specific hashtags<br><br>
<b>For CMO (CellPlex) demultiplexing - use non-AntibodyCapture files:</b><br>
• <b>10x_CMO</b>: CellPlex lipid barcodes (CMO301-CMO312)<br><br>
<b>For OCM (On-Chip Multiplexing):</b><br>
• Leave empty - no reference file needed<br><br>
<i>Contact genomics team to add custom barcode files.</i>"
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
    @params['customProbesFile', 'description'] = 'Custom probeset CSV-file according to 10X specifications (https://tinyurl.com/10xProbeSetCSVFormat). ONLY for probe-based single cell fixed RNA profiling (FRP). Note that all genes listed must have a corresponding entry in secondRef or controlSeqs. Custom probes must have the same length as the probes in the reference file.'
    @params['FeatureBarcodeFile'] = ''
    @params['FeatureBarcodeFile', 'file_upload'] = true
    @params['FeatureBarcodeFile', 'description'] = "Upload feature reference CSV for CITE-seq/ADT experiments.<br>
Required columns: <code>id, name, read, pattern, sequence, feature_type</code><br>
<a href='https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-feature-ref-csv'>10x Format Guide</a>"
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
