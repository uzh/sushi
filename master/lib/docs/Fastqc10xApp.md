# Fastqc10xApp

Read-quality control for 10x Genomics single-cell libraries.
ezRun backend: `EzAppFastqc_10x` (`R/app-fastQC_10x.R`).

## Scope

- **Analysis type:** quality control of the sequencing reads of 10x single-cell
  libraries, before CellRanger or any other processing.
- **Data:** 10x libraries delivered as one tar archive per sample (the `RawDataDir`
  column), always read as paired-end (R1 and R2).
- **Output:** quality-control reports only. The reads themselves go to later apps as
  they were.
- **Not covered here:** bulk sequencing reads delivered as FASTQ files
  (`FastqcApp`); screening reads against other genomes (`FastqScreen10xApp`).
- **Place in a pipeline:** a final QC step. No other app takes its output as input.

## What it can do

- **FastQC on one read pair per library.** From each sample's tar archive, only the
  first R1 file and the first R2 file are extracted and assessed. A library whose
  archive holds several FASTQ files per read (for example one per lane) is therefore
  assessed on part of its reads, not all of them.
- **FastQC's checks:** per-base and per-sequence quality, per-tile quality, per-base
  sequence content, GC content, N content, sequence length distribution,
  duplication, overrepresented sequences and adapter content. Adapters are searched
  with an in-house FGCZ adapter list.
- **Native report.** An overview page per FastQC plot type, a summary table of the
  FastQC checks per file, per-base quality heatmaps computed with ShortRead from a
  sample of reads per file, and the read counts from the dataset.
- **MultiQC aggregation.** All FastQC results are also combined into one MultiQC
  report.

## How to use it

- **Input dataset:** needs `Name`, `RawDataDir` (a `.tar` per sample) and
  `Read Count`. The run stops if `Read Count` is missing or an input is not a tar.
- **Jobs:** one job for the whole dataset.
- **Parameters:**

  | Parameter | Effect |
  |---|---|
  | `paired` | Shown on the form, but the app always treats the data as paired-end. |
  | `cmdOptions` | Extra options passed as-is to FastQC. |
  | `cores`, `ram`, `scratch` | Compute resources. |

- **Outputs:**
  - `Html [Link]`: the native summary report (`00index.html`).
  - `MultiQC [Link]`: the MultiQC report (`multiqc_report.html`).
  - `Report [File]`: the report folder, with the per-file FastQC results.
- **Gotchas:**
  - Only the first R1/R2 pair of each archive is assessed (see above).
  - Two extracted files with the same base name stop the run.

## What to describe in the Methods

A Methods description of this step is short, usually two or three sentences.

**Worth describing**

- **FastQC and its version.** FastQC does not print its version; it is the one in the
  job script's `module add` line (`QC/FastQC/<version>`).
- **Which reads were assessed:** the first R1/R2 file pair of each library, as
  paired-end reads. The files on the FastQC command line show which ones.
- **What FastQC assessed:** the checks listed under "What it can do".
- **The in-house adapter list** used for adapter content. Name it as such, not by
  its file path.
- **Options passed through `cmdOptions`,** where they change what FastQC assesses.
- **MultiQC and its version,** used to aggregate the results. MultiQC runs without
  additional options. The version is on its banner line
  (`/// MultiQC ... v<version>`). A `version_check ... now available!` line only
  names a newer release, not the one that ran.

**Not worth describing**

- **FastQC's housekeeping options:** `--extract`, `-o`, `-t` (threads) and the
  redirection of its output. They do not change the result.
- **K-mers.** The `--kmers` option has no effect, because FastQC's k-mer module is
  switched off in the FGCZ installation.
- **The native report and its ShortRead quality heatmaps:** presentation of the same
  QC.
- **Picard and samtools:** loaded, but never called on FASTQ input.
- **Tar extraction and file clean-up.**
- **Counts** of files, reports or reads. These are results.
- **R, ezRun and SUSHI.** They run the tools; they do not analyse the data.

**Example Methods paragraph** (base case, placeholders in angle brackets):

> Read quality of the 10x libraries was assessed with FastQC <version> on the first
> paired-end read files of each library, covering per-base and per-sequence quality,
> GC and N content, sequence length, duplication, overrepresented sequences and
> adapter content, with an in-house list of adapter sequences. The results were
> summarised with MultiQC <version>.
