# CellRangerMultiApp

Primary processing of 10x Genomics single-cell libraries with `cellranger multi`.
ezRun backend: `EzAppCellRangerMulti` (`R/app-cellRangerMulti.R`), with reference
helpers in `R/app-cellRanger.R`.

## Scope

- **Analysis type:** processing of 10x single-cell libraries from FASTQ to per-sample
  count matrices, including immune repertoire, antibody capture and sample
  demultiplexing where those libraries are present.
- **Designs covered:**
  - Gene expression: 3' and 5'.
  - T-cell receptor (VDJ-T) and B-cell receptor (VDJ-B).
  - Antibody capture: CITE-seq / ADT (Feature Barcoding).
  - Sample multiplexing: on-chip (OCM), hashtag antibodies (HTO), CellPlex (CMO).
  - Flex (fixed RNA, probe-based), including multiplexed Flex.
- **Not covered here:**
  - ATAC and multiome (`CellRangerATACApp`, `CellRangerARCApp`).
  - Spatial data (`SpaceRangerApp`).
  - Downstream analysis beyond CellRanger's own (`ScSeuratApp`).
  - Plain gene expression can also run through `CellRangerApp` (`cellranger count`).
- **Place in a pipeline:** after demultiplexing. Its per-sample output rows are the
  input to `CellBenderApp` and `ScSeuratApp`.

## What it can do

- **Any combination of the library types above in one run.** Each library comes as a
  tar archive or as FASTQ files. Several tars from one chip can be combined into one
  run.
- **Reference handling:**
  - By default, a shared CellRanger index for the selected genome, annotation and
    transcript types. It is reused when it exists, and built on first use with
    `cellranger mkref` after reducing the annotation to the chosen transcript types.
  - With `secondRef` or `controlSeqs`, a per-run custom reference: the genome plus
    extra sequences, for example transgenes or viral genes. Gene 3' ends can also be
    extended.
  - VDJ runs use a shared VDJ index, built on first use with `cellranger mkvdjref`.
- **Flex:**
  - The selected 10x probe set is reduced to probes whose genes are in the reference
    annotation, and custom probes can be added.
  - The chemistry is inferred when not given.
- **Multiplexing:**
  - HTO and CMO barcode sets are reduced to the barcodes in use.
  - Antibody capture and HTO share one feature reference.
  - Sample-to-barcode assignments come from Sample2Barcode files.
- **Inside `cellranger multi`:** read alignment and counting, cell calling, and,
  depending on the library types, antibody counting, V(D)J assembly and clonotype
  grouping, and assignment of cells to samples. Whenever gene expression is present,
  CellRanger's secondary analysis also runs (PCA, UMAP, clustering, differential
  expression between clusters); this app has no option to switch it off.
- **BAM handling:** BAM files are kept or deleted (`keepBam`). When kept and
  `secondRef` is set, they are stored as CRAM.
- **Read subsampling (hidden option).** `nReads` through `specialOptions`, for tar
  input only.

## How to use it

- **Input dataset.** One row per library or pool; one job per row.

  | Library type | Tar column | FASTQ columns |
  |---|---|---|
  | GEX / Flex | `RawDataDir` | `Read1`, `Read2` |
  | VDJ-T | `VdjTDataDir` | `VdjT Read1`, `VdjT Read2` |
  | VDJ-B | `VdjBDataDir` | `VdjB Read1`, `VdjB Read2` |
  | Feature Barcoding / HTO | `FeatureDataDir` | `Feature Read1`, `Feature Read2` |
  | CMO | `MultiDataDir` | `Multi Read1`, `Multi Read2` |

- **Sample `Name`:** must equal the GEX FASTQ file prefix. Otherwise CellRanger's
  preflight fails within a minute.
- **Multiplexed runs:** also need `<prefix>_Sample2Barcode.csv` files in the order's
  `o<order>_metaData/` folder on gStore. The prefix must match the start of the
  sample `Name`. The file can be generated with the Sample2Barcode Shiny app.
- **Dataset curation:** anything beyond plain gene expression usually needs the input
  dataset curated by hand.
- **Parameters:**

  | Parameter | Effect |
  |---|---|
  | `TenXLibrary` | Library types to process (several allowed). |
  | `MultiplexingType` | OCM, HTO (`antibody`) or CMO. Not used for Flex. |
  | `MultiplexBarcodeSet` | Barcode reference for HTO or CMO. |
  | `probesetFile`, `customProbesFile` | Flex probe set, and added custom probes. |
  | `FeatureBarcodeFile` | Antibody (ADT) feature reference. Not for HTO. |
  | `refBuild`, `refFeatureFile`, `transcriptTypes` | Reference genome, annotation and transcript types. |
  | `secondRef`, `controlSeqs` | Extra sequences, which trigger a per-run custom reference. |
  | `includeIntrons` | Count intronic reads. Ignored for Flex. |
  | `chemistry` | Pin the chemistry instead of auto-detection. |
  | `expectedCells` | Expected cell number. Empty lets CellRanger estimate it. |
  | `keepBam` | Whether the alignment files are kept. |
  | `cmdOptions` | Extra options passed as-is to `cellranger multi`. |
  | `CellRangerVersion` | CellRanger version. It also changes the output folder layout. |
  | `cores`, `ram`, `scratch` | Compute resources. Large multiplexed Flex runs need more memory. |

- **Outputs:**
  - `ResultDir [File,Link]`: CellRanger's `outs` folder, renamed to the sample name.
    It contains `config.csv` and, per demultiplexed sample,
    `per_sample_outs/<sample>/`.
  - `Report [Link]`: the per-sample web summary (CellRanger 9 and earlier) or the
    combined QC report (10 and later).
  - One further row per demultiplexed sample, **only for multiplexed runs**, carrying
    `CountMatrix`, `UnfilteredCountMatrix`, `ResultDir` and `Report`.
- **Gotchas:**
  - The main output row has no `CountMatrix`. A run without multiplexing therefore
    produces no row that `ScSeuratApp` or `CellBenderApp` can take directly.
  - `ScSeuratApp` also needs a `Condition` column, present only when the input
    dataset carries one.
  - Dropdown parameters left at `select` break submissions made outside the web form.
  - A run waits while another job builds the same shared reference, and gives up after
    a fixed timeout.
  - Operating and debugging notes are in the `cellranger-fgcz` and
    `sample2barcode-generation` skills.

## What to describe in the Methods

Most of the analysis happens inside `cellranger multi`. The `config.csv` in the
result folder is the configuration CellRanger actually received. Where it and the
job script differ, `config.csv` is what applied.

**Which internal steps run.** They depend on the library types, as listed in
`config.csv`:

| In `config.csv` | Steps performed inside CellRanger |
|---|---|
| `Gene Expression` in `[libraries]`, no `probe-set` line | Alignment of reads to the reference genome, UMI counting per gene, cell calling; secondary analysis (PCA, UMAP, clustering, differential expression between clusters) |
| `Gene Expression` in `[libraries]`, with a `probe-set` line (Flex) | Alignment of reads to the probe set, UMI counting per gene, cell calling; secondary analysis as above |
| `VDJ-T` or `VDJ-B` in `[libraries]` | Assembly and annotation of T-cell or B-cell receptor contigs, and grouping into clonotypes |
| `Antibody Capture` in `[libraries]`, and the `[feature]` reference is not a `multi_barcode_set...csv` file (that one holds only hashtags) | Counting of antibody-derived tags per cell |
| `[samples]` with `cmo_ids` | Assignment of cells to samples by CellPlex (CMO) tags |
| `[samples]` with `hashtag_ids` | Assignment of cells to samples by hashtag antibodies |
| `[samples]` with `ocm_barcode_ids` | Assignment of cells to samples by on-chip multiplexing barcodes |
| `[samples]` with `probe_barcode_ids` | Assignment of cells to samples by Flex probe barcodes |

CellRanger chooses the methods and settings of these steps internally.

**Worth describing**

- **Cell Ranger and its version.** The version is the one in the job script's
  `module load Aligner/CellRanger/<version>` line.
- **The library types processed,** from the `feature_types` in `[libraries]`.
- **The reference.** The gene expression reference path in `config.csv` follows
  `<organism>/<source>/<genome build>/Annotation/Release<_><release>-<date>/Genes/<index>`:
  - It gives the organism, the annotation source (e.g. Ensembl, GENCODE), the genome
    build and the annotation release.
  - The `-<date>` is when FGCZ set the reference up; it is not part of the release.
  - The index folder name `genes_10XGEX_SC_<types>_Index` lists the transcript types
    the annotation was reduced to.
  - A `10X_customised_Ref` folder means the genome was extended with extra sequences.
  - The VDJ reference follows the same layout.
- **Settings in `config.csv`,** when present: `chemistry`, `include-introns`,
  `expect-cells`, and the Flex probe set, whose file name gives the 10x probe set and
  version. For Flex, the probe set was restricted to genes present in the reference
  annotation.
- **Antibody panel and multiplexing:** that a user-supplied antibody panel was used,
  and the number of multiplexed samples listed in `[samples]`.
- **The internal steps from the table above,** for the matching rows only.
- **Reference building, when it happened.** A `cellranger mkref` or
  `cellranger mkvdjref` line means the reference was built during the run.
- **Read subsampling, when it happened:** `seqtk sample` lines, with the seed and
  read number on them.

**Not worth describing**

- **Housekeeping options of `cellranger multi`:** `--id`, `--localmem`,
  `--localcores` and `--csv`. Also `create-bam` and the scheduler resources.
- **Cell-type annotation.** CellRanger checks for it and skips it, because FGCZ
  references are not among those it supports for annotation; its log reports the
  skip. No cell types are produced.
- **`includeIntrons` in the job script of a Flex run:** it is not passed to
  CellRanger.
- **File handling:** tar extraction, FASTQ linking, renames, BAM deletion or CRAM
  conversion, the `expanded_dataset.tsv` and the copies to gStore.
- **seqtk and samtools:** loaded, but run only for subsampling and CRAM conversion.
- **Counts** of cells, reads or genes. These are results.
- **R, ezRun and SUSHI.** They run CellRanger; they do not analyse the data.

**Example Methods paragraph** (base case for a gene expression library, placeholders
in angle brackets; one sentence in the same style is added for each further matching
row of the table):

> Single-cell gene expression libraries were processed with Cell Ranger <version>
> (cellranger multi) against the <organism> <genome build> genome with <source>
> annotation release <release>, restricted to <transcript types>. Read alignment, UMI
> counting and cell calling were performed within Cell Ranger, which also ran its
> secondary analysis (principal component analysis, UMAP, clustering and
> differential expression between clusters).

A further sentence in the same style, for a `VDJ-T` library:

> T-cell receptor contigs were assembled, annotated and grouped into clonotypes within
> Cell Ranger.
