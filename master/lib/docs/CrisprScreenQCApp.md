# CrisprScreenQCApp

Quality control of CRISPR screen libraries: identification of the sgRNA library and
read composition.
ezRun backend: `EzAppCrisprScreenQC` (`R/app-CrisprScreenQCApp.R`), report template
`inst/templates/CrisprScreenQC.Rmd`.

## Scope

- **Analysis type:** a QC screen of CRISPR screen samples. A subsample of reads is
  counted against every sgRNA library in the FGCZ collection at once. The report shows
  which library each sample matches, the top sgRNAs, and the base composition of the
  reads.
- **Data:** FASTQ reads of pooled CRISPR screens, single-end or paired-end (only the
  first read is counted).
- **Output:** a report only.
- **Not covered here:** full counting against one library (`CountSpacerApp`,
  `MageckCountApp`); testing (`MageckTestApp`).
- **Place in a pipeline:** a QC step before counting. No other app takes its output
  as input.

## What it can do

- **Merged library reference.** The MAGeCK library files of all sgRNA libraries in the
  FGCZ collection (`libPath`) are merged into one reference, with each sgRNA labelled
  by its library.
- **Subsampling.** Each sample is reduced to `nReads` reads with ShortRead, using a
  seed fixed in the code.
- **Preprocessing with fastp.** Adapter trimming is always on, with the usual optional
  trimming and filtering.
- **Counting with `mageck count`** of each sample against the merged reference.
  MAGeCK determines the 5' trim length of the reads itself.
- **Library identification:** the library with the most counted reads is reported as
  the sample's match, with the top `topFeatures` sgRNAs.
- **Base composition:** per-position base frequencies of the trimmed reads, as
  sequence logos.

## How to use it

- **Input dataset:** needs `Name` and `Read1`; `paired` is set from the presence of
  `Read2`.
- **Jobs:** one job for the whole dataset.
- **Parameters:**

  | Parameter | Effect |
  |---|---|
  | `nReads` | Reads per sample used for the screen. |
  | `libPath` | The sgRNA library collection searched. |
  | `topFeatures` | Number of top sgRNAs shown per sample. |
  | fastp parameters | Read preprocessing. A value of 0 switches an option off. |
  | `cores`, `ram`, `scratch` | Compute resources. |

- **Outputs:** `Html [Link]` (the report) and `Report [File]`.

## What to describe in the Methods

A Methods description of this step is usually one or two sentences. It is a QC step.
MAGeCK prints its version ("Welcome to MAGeCK v<version>") and its full command
("Parameters: ... mageck count ...") in the log. fastp's command line is logged as
`EXECUTED CMD: fastp ...`.

**Worth describing**

- **The subsample:** `nReads` reads per sample.
- **fastp preprocessing,** with its version and the options it was given that are
  switched on, including adapter trimming. fastp also applies its built-in quality
  filter, and may trim poly-G tails on data from two-colour Illumina instruments.
- **MAGeCK and its version,** used to count the reads against the combined sgRNA
  libraries of the FGCZ collection, with automatic detection of the 5' trim length.
- **Library identification** by the highest read count.

**Not worth describing**

- **The merging of the library files and the report layout.**
- **The base composition logos:** presentation only.
- **fastp options set to 0,** which are switched off, and fastp's `--thread` and
  `--compression`.
- **Counts** of reads, sgRNAs or matches. These are results.
- **R, ezRun and SUSHI.** They run the tools; they do not analyse the data.

**Example Methods paragraph** (placeholders in angle brackets):

> For each sample, <nReads> reads were subsampled, preprocessed with fastp <version>
> (<switched-on options>) and counted with MAGeCK <version> (mageck count) against the
> combined sgRNA libraries of the FGCZ collection, to identify the library used.
