# CellRangerApp

Primary processing of 10x Genomics single-cell libraries with `cellranger count`
(gene expression) or `cellranger vdj` (immune repertoire).
ezRun backend: `EzAppCellRanger` (`R/app-cellRanger.R`).

## Scope

- **Analysis type:** processing of one 10x library per job, from FASTQ to a
  cell-by-gene count matrix (gene expression), or to assembled T- or B-cell receptor
  contigs (VDJ).
- **Data:** 10x single-cell libraries without multiplexing, antibody capture or Flex:
  3' or 5' gene expression, or a VDJ library alone.
- **Output:** CellRanger's output folder per sample, with count matrices, a web
  summary, and optionally the alignments as CRAM.
- **Not covered here:**
  - multiplexed, antibody-capture, Flex, or combined GEX+VDJ designs
    (`CellRangerMultiApp`)
  - ATAC and multiome (`CellRangerATACApp`, `CellRangerARCApp`)
  - spatial data (`SpaceRangerApp`)
- **Place in a pipeline:** after demultiplexing. Its output rows feed `CellBenderApp`
  and `ScSeuratApp`.

## What it can do

- **Input staging.** FASTQ files are linked, or tar archives extracted. Several tar
  archives of one sample (for example two sequencing runs) are combined. Tar mode
  only: read subsampling with seqtk (`nReads` through `specialOptions`).
- **Gene expression reference:**
  - By default, a shared CellRanger index for the selected genome, annotation and
    transcript types. It is reused when it exists, and built on first use with
    `cellranger mkref` after reducing the annotation to the chosen transcript types.
  - With `secondRef` or `controlSeqs`, a per-run custom reference: the genome plus
    extra sequences.
- **`cellranger count`** (`TenXLibrary = GEX`):
  - alignment, UMI counting per gene and cell calling, with the chemistry set or
    auto-detected, optional intronic counting (`includeIntrons`), and an optional
    expected cell number
  - CellRanger's secondary analysis always runs: PCA, UMAP, clustering and
    differential expression between clusters
- **Cell-type annotation (human, CellRanger 10.1 and later).** CellRanger's local
  Azimuth pan-human model annotates cells when the reference is recognised as human.
  For human references, ezRun passes CellRanger an alias of the reference that
  declares the genome as GRCh38, so the annotation runs; the results are in
  `cell_types/`. For other organisms, or older CellRanger versions, CellRanger skips
  annotation and says so in its log.
- **`cellranger vdj`** (`TenXLibrary = VDJ`): contig assembly, annotation and
  clonotype grouping, against a shared VDJ index built with `cellranger mkvdjref` on
  first use.
- **Optional extras:**
  - `runVeloCyto`: spliced and unspliced counts with velocyto, for RNA velocity
  - `bamStats`: per-cell alignment statistics
  - `keepAlignment`: keeps the alignments, converted to CRAM when the reference is
    the standard one
- **Output layout.** CellRanger's `outs` folder is renamed to the sample name.

## How to use it

- **Input dataset:** needs `Name` and `Species`, plus either `RawDataDir` (tar
  archives) or `Read1` and `Read2`. The sample `Name` must match the FASTQ file prefix.
- **Jobs:** one job per sample.
- **Parameters:**

  | Parameter | Effect |
  |---|---|
  | `TenXLibrary` | `GEX` (`cellranger count`) or `VDJ` (`cellranger vdj`). |
  | `refBuild`, `refFeatureFile`, `transcriptTypes` | Reference genome, annotation, and the transcript types kept. |
  | `chemistry` | Chemistry, or `auto` for detection by CellRanger. |
  | `includeIntrons` | Count intronic reads (default true). |
  | `expectedCells` | Expected cell number; empty lets CellRanger estimate it. |
  | `controlSeqs`, `secondRef` | Extra sequences, which trigger a per-run custom reference. |
  | `runVeloCyto`, `bamStats`, `keepAlignment` | Optional extras, as above. |
  | `cmdOptions` | Extra options passed as-is to CellRanger. |
  | `CellRangerVersion` | CellRanger version. |
  | `cores`, `ram`, `scratch` | Compute resources. |

- **Outputs (per sample):**
  - `ResultDir [File]`: the CellRanger output folder.
  - `Report [Link]`: the web summary.
  - `CountMatrix [Link]`, `UnfilteredCountMatrix [Link]`: filtered and raw matrices
    (GEX).
  - `AlignmentFile [Link]`: the CRAM file, when kept.

## What to describe in the Methods

The CellRanger command line is logged as `EXECUTED CMD: cellranger count ...` (or
`vdj`), and the CellRanger version is the job script's
`module load Aligner/CellRanger/<version>` line.

**Worth describing**

- **Cell Ranger and its version,** with `count` or `vdj`.
- **The reference.** The `--transcriptome` (or `--reference`) path follows
  `<organism>/<source>/<genome build>/Annotation/Release<_><release>-<date>/Genes/<index>`:
  - It gives the organism, the annotation source (e.g. Ensembl, GENCODE), the genome
    build and the annotation release.
  - The `-<date>` is when FGCZ set the reference up; it is not part of the release.
  - The index folder name `genes_10XGEX_SC_<types>_Index` lists the transcript types
    the annotation was reduced to.
  - A `10X_annotatable_Ref` folder is the same reference under an alias, used so
    that CellRanger annotates human cells.
  - A `10X_customised_Ref` folder means the genome was extended with extra sequences.
  - A `cellranger mkref` or `mkvdjref` line means the reference was built during the
    run.
- **The settings on the command line, in words:** the chemistry (`auto` means
  detected by CellRanger), intronic reads included or not (`--include-introns`), the
  expected cell number when given, and options from `cmdOptions`.
- **The steps performed inside CellRanger:**
  - for `count`: alignment to the genome, UMI counting per gene and cell calling,
    followed by its secondary analysis (PCA, UMAP, clustering and differential
    expression between clusters)
  - for `vdj`: contig assembly and annotation, and grouping into clonotypes

  CellRanger chooses the methods and settings of these steps internally.
- **Cell-type annotation, when it ran:** a `cell_types/` folder in the output, with
  CellRanger's local Azimuth pan-human model.
- **RNA velocity counting with velocyto, when `runVeloCyto` was true.**
- **Read subsampling, when it happened:** `seqtk sample` lines.

**Not worth describing**

- **Housekeeping options:** `--id`, `--fastqs`, `--sample`, `--localmem`,
  `--localcores` and `--create-bam`.
- **A skipped cell-type annotation.** The log line "Cell annotation via 10x Cloud was
  skipped because the reference used is not yet supported" means no annotation was
  done.
- **File handling:** tar extraction, FASTQ linking, renames, BAM-to-CRAM conversion
  and deletions.
- **seqtk and samtools,** unless their own commands appear.
- **Counts** of cells, reads or genes. These are results.
- **R, ezRun and SUSHI.** They run CellRanger; they do not analyse the data.

**Example Methods paragraph** (gene expression base case, placeholders in angle
brackets):

> Single-cell gene expression libraries were processed with Cell Ranger <version>
> (cellranger count) against the <organism> <genome build> genome with <source>
> annotation release <release>, restricted to <transcript types>, with the chemistry
> <detected automatically | set to ...> and intronic reads <included | excluded>.
> Read alignment, UMI counting and cell calling were performed within Cell Ranger,
> which also ran its secondary analysis (principal component analysis, UMAP,
> clustering and differential expression between clusters).
