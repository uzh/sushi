# SpaceRangerApp

Primary processing of 10x Genomics Visium spatial data with `spaceranger count`.
ezRun backend: `EzAppSpaceRanger` (`R/app-spaceRanger.R`), with reference helpers in
`R/app-cellRanger.R`.

## Scope

- **Analysis type:** processing of one Visium capture area per job, from FASTQ and
  tissue image to spatially resolved count matrices. The tissue is detected and the
  reads are aligned to the image by fiducials or a manual alignment. For Visium HD,
  counts come in several bin sizes, with optional cell segmentation.
- **Data:** Visium (probe-based or poly-A), Visium CytAssist and Visium HD, for human
  or mouse probe sets or any organism with an FGCZ reference.
- **Output:** SpaceRanger's output folder per sample, with count matrices, the
  spatial images, the web summary and a pseudo-bulk count table.
- **Not covered here:** Xenium (`XeniumQCApp`, `XeniumSeuratApp`); downstream
  spatial analysis (`SpatialSeuratApp`, `SpatialSeuratHDApp`, `VisiumHDSeuratApp`).
- **Place in a pipeline:** after demultiplexing. Its outputs feed the spatial Seurat
  apps and `VisiumQCApp`.

## What it can do

- **Input staging.** FASTQ files are linked, or tar archives extracted; several
  archives of one sample are combined.
- **Reference:**
  - By default, a shared SpaceRanger/CellRanger index for the selected genome,
    annotation and transcript types, built on first use.
  - With `secondRef` or `controlSeqs`, a per-run custom reference.
- **Probe sets (`probesetFile`).** The 10x probe set is reduced to probes whose genes
  are in the reference annotation, and custom probes can be added
  (`customProbesFile`).
- **Images:**
  - the brightfield image (`Image`), a dark or fluorescence image (`darkImage`), and
    the CytAssist image (`CytaImage`)
  - for Visium HD with a multi-page TIFF under 4 GB, the highest-resolution page is
    extracted and used
  - a manual Loupe alignment file (`loupe-alignment`) replaces automatic fiducial
    alignment
- **Slide and area.** Taken from the dataset (`Slide`, `Area`), so SpaceRanger uses
  the slide's layout file. Without them, SpaceRanger runs in unknown-slide mode.
- **Antibody or gene panels (`panelFile`).** An additional panel library, with its
  feature reference.
- **Inside `spaceranger count`:**
  - tissue detection and image alignment
  - read alignment (to the genome, or to the probe set for probe-based assays) and
    UMI counting per spot or bin
  - secondary analysis: PCA, UMAP, clustering and differential expression between
    clusters. For Visium HD this runs at 8 and 16 um bins, next to the 2 um
    matrices.
- **Visium HD cell segmentation (`runSegmentation`).** Nucleus and cell segmentation
  from the image, giving per-cell matrices and their own secondary analysis. When
  off, segmentation is switched off.
- **Cell annotation.** SpaceRanger's own cell annotation can be switched off through
  `cmdOptions` (`--disable-cell-annotation`).
- **After the run:**
  - the alignments are kept as CRAM or deleted (`keepAlignment`)
  - a pseudo-bulk table is written: gene counts summed over all filtered spots, or
    over the 16 um bins for Visium HD

## How to use it

- **Input dataset:** needs `Name`, `Species`, `Slide` and `Area`, plus either
  `RawDataDir` (tar archives) or `Read1` and `Read2`. Optional columns: `Image`,
  `CytaImage`, `loupe-alignment`, and `PanelRawDataDir` for a panel library.
- **Jobs:** one job per capture area.
- **Parameters:**

  | Parameter | Effect |
  |---|---|
  | `refBuild`, `refFeatureFile`, `transcriptTypes` | Reference genome, annotation, and the transcript types kept. |
  | `probesetFile`, `customProbesFile` | 10x probe set for probe-based assays, and added custom probes. |
  | `panelFile` | Feature reference of an additional panel library. |
  | `controlSeqs`, `secondRef` | Extra sequences, which trigger a per-run custom reference. |
  | `runSegmentation` | Visium HD nucleus and cell segmentation. |
  | `keepAlignment` | Keep the alignments as CRAM. |
  | `cmdOptions` | Extra options passed as-is to SpaceRanger. |
  | `SpaceRangerVersion` | SpaceRanger version. |
  | `cores`, `ram`, `scratch` | Compute resources. |

- **Outputs (per sample):** the SpaceRanger output folder (renamed to the sample),
  with `filtered_feature_bc_matrix`, `spatial/`, `web_summary.html`, for Visium HD
  `binned_outputs/` and `segmented_outputs/`, and the pseudo-bulk `-counts.txt`.

## What to describe in the Methods

The SpaceRanger command line is logged as `EXECUTED CMD: spaceranger count ...`, and
the SpaceRanger version is the job script's `module load Aligner/SpaceRanger/<version>`
line.

**Worth describing**

- **Space Ranger and its version.**
- **The assay and slide:** Visium, CytAssist or Visium HD, with the slide and capture
  area (`--slide`, `--area`), or unknown-slide mode.
- **The images used:** brightfield (`--image`), dark or fluorescence
  (`--darkimage`), CytAssist (`--cytaimage`), and a manual alignment
  (`--loupe-alignment`) when given.
- **The reference.** The `--transcriptome` path follows
  `<organism>/<source>/<genome build>/Annotation/Release<_><release>-<date>/Genes/<index>`:
  - It gives the organism, the annotation source, the genome build and the
    annotation release.
  - The `-<date>` is when FGCZ set the reference up; it is not part of the release.
  - The index folder name `genes_10XGEX_SC_<types>_Index` lists the transcript types
    included.
- **The probe set,** from the `--probe-set` file name (10x probe set and version),
  restricted to genes present in the reference annotation, with custom probes when
  given.
- **A panel library, when `--feature-ref` is present.**
- **The steps performed inside Space Ranger:**
  - tissue detection and image alignment
  - read alignment and UMI counting per spot or bin
  - its secondary analysis (PCA, UMAP, clustering and differential expression between
    clusters)
  - for Visium HD: counting at 2, 8 and 16 um bins, and nucleus and cell segmentation
    when it ran (a `segmented_outputs/` folder)

  Space Ranger chooses the methods and settings of these steps internally.
- **Options from `cmdOptions`** that change the analysis, for example
  `--disable-cell-annotation`.

**Not worth describing**

- **Housekeeping options:** `--id`, `--fastqs`, `--sample`, `--localmem`,
  `--localcores` and `--create-bam`.
- **The TIFF page extraction** (`tiffsplit`), which only picks the high-resolution
  image.
- **File handling:** tar extraction, FASTQ linking, renames, BAM-to-CRAM conversion
  and deletions.
- **The pseudo-bulk table,** a summary for convenience.
- **The CellRanger and samtools modules:** support only.
- **Counts** of spots, bins, cells, reads or genes. These are results.
- **R, ezRun and SUSHI.** They run Space Ranger; they do not analyse the data.

**Example Methods paragraph** (Visium HD base case, placeholders in angle brackets):

> Visium HD data were processed with Space Ranger <version> (spaceranger count) using
> the <probe set name and version>, restricted to genes of the <source> annotation
> release <release> for <organism> <genome build>, with the high-resolution
> brightfield image and the CytAssist image, for slide <slide>, area <area>. Tissue
> detection, image alignment, read alignment to the probe set and UMI counting were
> performed within Space Ranger, which also ran its secondary analysis (principal
> component analysis, UMAP, clustering and differential expression) at 8 and 16 um
> bin sizes and segmented nuclei and cells from the image.
