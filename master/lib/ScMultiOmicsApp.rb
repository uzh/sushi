#!/usr/bin/env ruby
# encoding: utf-8

# ScMultiOmicsApp
# ----------------
# Downstream multi-omics extension layered on top of an annotated ScSeurat
# scData.qs2 (or a BD Rhapsody pre-built Seurat). For one sample per job
# (process_mode = SAMPLE), the runtime auto-detects which modalities are
# present alongside the original CountMatrix and produces a single multi-assay
# Seurat object plus an HTML report.
#
# Detection rules (R side, see ezRun::detectModalities):
#   - RNA   : filtered_feature_bc_matrix.h5 next to CountMatrix
#   - ADT   : H5 carries any "Antibody Capture" feature_type
#   - VDJ-T : sibling vdj_t/filtered_contig_annotations.csv (CellRanger Multi)
#   - VDJ-B : sibling vdj_b/filtered_contig_annotations.csv (CellRanger Multi)
#   - ATAC  : sibling atac_fragments.tsv.gz + atac_peaks.bed (CellRanger ARC)
#
# Output dataset surfaces:
#   - Report [File]      : the per-sample report directory on gstore
#   - Static Report [Link] : 00index.html
#   - ScMultiOmics [Link]  : scMultiData.qs2 (Seurat with DefaultAssay = RNA)
#
# DefaultAssay is kept on RNA so the existing exploreSC Shiny app keeps
# rendering the RNA UMAP unchanged.

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class ScMultiOmicsApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'ScMultiOmics'
    # Single-sample app: one job per row in the input dataset.
    @params['process_mode'] = 'SAMPLE'
    @analysis_category = 'SingleCell'

    @description =<<-EOS
Downstream multi-omics extension on top of an annotated ScSeurat scData.qs2.<br/>
Auto-discovers ADT, VDJ, and ATAC siblings of the original CountMatrix and
produces a multi-assay Seurat object plus an HTML report. RNA assay stays
default so existing exploreSC apps render the RNA UMAP unchanged.<br/>
Supports CellRanger Multi (RNA + ADT + VDJ-T/B), CellRanger ARC
(RNA + ATAC), standalone VDJ via VDJTPath/VDJBPath columns, and BD Rhapsody
(SCDataOrigin = "BDRhapsody").
    EOS

    # Required dataset.tsv columns:
    #   Name       - sample name (used to prefix barcodes <Sample>_<barcode>)
    #   Species    - human-readable species (informational)
    #   refBuild   - genome reference path (used to resolve EnsDb annotation
    #                for Signac GeneActivity on ATAC; auto-detected human/mouse)
    #   CountMatrix - filtered matrix directory or H5 parent. Phase 1 expects
    #                the CellRanger filtered output (mtx dir or count/sample_*).
    #
    # Optional columns (auto-honored when present):
    #   Report / SC Cluster Report / SC Seurat
    #     - link to a previous ScSeurat job; we reuse its scData.qs2.
    #   VDJTPath / VDJBPath
    #     - directory containing filtered_contig_annotations.csv;
    #       overrides auto-discovery when set.
    #   BDRhapsodyPath
    #     - directory containing the BD Rhapsody *_Seurat.rds; together with
    #       SCDataOrigin = "BDRhapsody" routes the runtime through the BD loader.
    @required_columns = ['Name', 'Species', 'refBuild', 'CountMatrix']
    @required_params = ['name']

    # ---------------------------------------------------------------------
    # SLURM resources
    # ---------------------------------------------------------------------
    # Cores feed Seurat's `future::plan("multicore")` and the qs2 nthreads
    # used to read scData.qs2 / write scMultiData.qs2. 4 cores is enough for
    # most 10x Multi or ARC samples; bump to 8+ for >100k cells.
    @params['cores'] = '4'
    @params['cores', "context"] = "slurm"
    # RAM in GB. ATAC fragments parsing with Signac is the main consumer;
    # 64 GB covers typical multiome samples up to ~30k cells. Increase for
    # ARC samples > 50k cells or BD Rhapsody runs (which can exceed 60k cells
    # and carry ADT panels with hundreds of markers).
    @params['ram'] = '64'
    @params['ram', "context"] = "slurm"
    # Local scratch in GB for staging intermediates (qs2 reads/writes,
    # ATAC fragment scans). 50 GB is sufficient for every dataset class
    # currently supported.
    @params['scratch'] = '50'
    @params['scratch', "context"] = "slurm"

    @params['name'] = 'ScMultiOmics'
    # Reference genome assembly. The R runtime parses this string to pick the
    # right EnsDb package for ATAC GeneActivity (Hsapiens for GRCh38/hg38,
    # Mmusculus for GRCm/mm10/mm39). When neither matches, GeneActivity is
    # skipped silently and the rest of the ATAC tab still renders.
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "reference genome assembly"
    # Carried through so downstream ScSeuratCombine (which gates on
    # refFeatureFile in @required_columns) can consume our output unchanged.
    @params['refFeatureFile'] = 'genes.gtf'
    @params['refFeatureFile', "context"] = "ScMultiOmics"

    # ---------------------------------------------------------------------
    # ADT options
    # ---------------------------------------------------------------------
    @params['adtNorm', 'hr-header'] = "ADT"
    # ADTnorm = peak-deconvolution normalization recommended for marker-rich
    # panels (TotalSeq B/C, BD AbSeq) - slow but more biologically meaningful.
    # CLR = Seurat default centered-log-ratio - fast, good for sanity checks.
    @params['adtNorm'] = ['ADTnorm', 'CLR']
    @params['adtNorm', 'description'] = "ADT normalization method. ADTnorm is recommended for marker-rich panels; CLR is faster and matches Seurat's default."
    # ADT PC count. Capped at (#features - 1) at runtime so small panels
    # don't fail the SVD step.
    @params['npcsADT'] = 18
    @params['npcsADT', 'description'] = "Number of PCs to compute on the ADT assay (clamped to #ADT features - 1 at runtime)."

    # ---------------------------------------------------------------------
    # VDJ options
    # ---------------------------------------------------------------------
    @params['vdjChain', 'hr-header'] = "VDJ"
    # 'auto' picks whichever chains the modality detector found. 'TCR' or
    # 'BCR' force a specific chain even when the other is also available.
    # 'both' currently expects sequential TCR + BCR loading; v1 only attaches
    # one assay so 'both' falls back to TCR.
    @params['vdjChain'] = ['auto', 'TCR', 'BCR', 'both']
    @params['vdjChain', 'description'] = "Which VDJ chains to load. 'auto' uses sibling vdj_t / vdj_b (CellRanger Multi) or VDJTPath / VDJBPath columns when set."

    # How TCR clones are merged inside scRepertoire::combineExpression.
    # 'strict' = V/J genes + CDR3 nucleotide (most stringent, recommended for
    #            within-donor clonal expansion).
    # 'aa'     = CDR3 amino acid sequence (merges synonymous mutations;
    #            useful for cross-donor convergence studies).
    @params['cloneCallTCR'] = ['strict', 'aa']
    @params['cloneCallTCR', 'description'] = "TCR clonotype definition: 'strict' (V/J + CDR3 nt, exact match) or 'aa' (CDR3 amino acid, merges synonymous mutations)."

    # Optional: collapse near-identical CDR3 sequences via
    # scRepertoire::clonalCluster (Hamming-style edit-distance clustering).
    # Off by default — most projects want exact-match clonality. Turn on for
    # cross-donor convergence studies or noise correction on long runs.
    @params['tcrSimilarityMerge'] = false
    @params['tcrSimilarityMerge', 'description'] = "Advanced: collapse near-identical TCR CDR3 sequences via scRepertoire::clonalCluster before counting clones. Off by default."

    @params['tcrSimilarityThreshold'] = 0.85
    @params['tcrSimilarityThreshold', 'description'] = "Normalized similarity (0-1) for TCR clonalCluster. Only effective when tcrSimilarityMerge = true. 0.85 ≈ 15% Hamming distance; tighten to 0.95 for noise-only correction."

    # BCR similarity merging is required for somatic-hypermutation lineages and
    # is on by default in scRepertoire (combineBCR threshold = 0.85). Surfaced
    # so users can tighten it on clean datasets or loosen it on heavily mutated
    # repertoires.
    @params['bcrSimilarityThreshold'] = 0.85
    @params['bcrSimilarityThreshold', 'description'] = "Similarity threshold for BCR clone merging (handles SHM). Tighten (0.90-0.95) for clean data; loosen (0.75) for highly mutated lineages."

    # ---------------------------------------------------------------------
    # WNN options
    # ---------------------------------------------------------------------
    @params['runWNN', 'hr-header'] = "WNN"
    # WNN runs whenever 2+ dimensional reductions exist on the object
    # (any subset of {pca, adt.pca, lsi}). VDJ is metadata-only and never
    # part of WNN. Disable to skip the integration step in single-modality
    # corner cases.
    @params['runWNN'] = true
    @params['runWNN', 'description'] = "Run Seurat::FindMultiModalNeighbors when 2+ dimensional modalities are present (RNA + ADT, RNA + ATAC, or RNA + ADT + ATAC)."

    # Resolution passed to Seurat::FindClusters on the WNN graph (algorithm 3,
    # Leiden). Higher = more, finer clusters; lower = fewer, coarser.
    @params['wnnResolution'] = 0.5
    @params['wnnResolution', 'description'] = "Resolution for FindClusters on the weighted SNN graph. Stored on the object as wsnn_res.<resolution>."

    # ---------------------------------------------------------------------
    # Misc
    # ---------------------------------------------------------------------
    @params['specialOptions'] = ''
    @params['mail'] = ''
    @modules = ["Dev/R"]
    # Inheriting Factor / B-Fabric tags propagates Condition columns and
    # B-Fabric metadata from the upstream dataset.
    @inherit_tags = ["Factor", "B-Fabric"]
  end

  def preprocess
    # Random suffix used by SUSHI to disambiguate result directories.
    @random_string = (1..12).map{[*('a'..'z')].sample}.join
  end

  def next_dataset
    # Output dataset row written to dataset.tsv after the job finishes.
    # The R runtime writes scMultiData.qs2 + 00index.html into report_file.
    report_file = File.join(@result_dir, "#{@dataset['Name']}_ScMultiOmicsReport")
    report_link = File.join(report_file, '00index.html')
    {'Name'=>@dataset['Name'],
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'refFeatureFile'=>@params['refFeatureFile'],
     'Static Report [Link]'=>report_link,
     'Report [File]'=>report_file,
     # ScMultiOmics [Link] is the path that exploreSC reads to render the
     # RNA UMAP (DefaultAssay = RNA in the saved object).
     'ScMultiOmics [Link]'=>File.join(report_file, "scMultiData.qs2"),
     # SC Seurat [Link] -> scData.qs2 (symlink to scMultiData.qs2 created in
     # ezMethodScMultiOmics). Lets ScSeuratCombine + ScSeuratCombinedLabelClusters
     # consume the multi-omics output as if it were a standard ScSeurat result.
     'SC Seurat [Link]'=>File.join(report_file, "scData.qs2"),
    }.merge(extract_columns(@inherit_tags))
  end

  def set_default_parameters
    # Inherit refBuild + refFeatureFile from the input dataset row by default.
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('refFeatureFile')
      @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    end
  end

  def commands
    # Defers to ezRun::EzAppScMultiOmics (R6 wrapper around ezMethodScMultiOmics).
    # See ~/git/ezRun/R/app-ScMultiOmics.R for the full runtime.
    run_RApp("EzAppScMultiOmics")
  end
end

if __FILE__ == $0

end
