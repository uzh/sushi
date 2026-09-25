# FastqcApp

Read-quality control for bulk sequencing data.
ezRun backend: `EzAppFastqc` (`R/app-fastQC.R`).

## Scope

- **Analysis type:** quality control of sequencing reads (FASTQ), before any alignment
  or quantification.
- **Data:** bulk (non-10x) short-read sequencing, single-end or paired-end, from any
  organism and any assay. No reference genome is involved.
- **Output:** quality-control reports only. The reads themselves go to later apps as
  they were.
- **Not covered here:** 10x single-cell reads (`Fastqc10xApp`); screening reads
  against other genomes for contamination (`FastqScreenApp`); trimming reads for
  downstream analysis.
- **Place in a pipeline:** a final QC step. No other app takes its output as input.

## What it can do

- **FastQC on every read file.** FastQC checks per-base and per-sequence quality,
  per-tile quality, per-base sequence content, GC content, N content, sequence
  length distribution, duplication, overrepresented sequences and adapter content.
  Adapters are searched with an in-house FGCZ adapter list.
- **MultiQC aggregation.** All FastQC results are combined into one MultiQC report.
- **Samples with several files.** A sample whose `Read1`/`Read2` cell lists several
  FASTQ files, comma-separated, has them concatenated first.
- **Subsampling of very large datasets.** When the total read count of the dataset
  is above a fixed threshold, each sample is reduced to a fixed number of reads
  before FastQC runs. The QC then describes the subsample, not all reads.
- **Native FastQC reports (optional).** With `showNativeReports`, the per-file FastQC
  reports are kept, and an overview page per plot type plus a summary report are
  rendered.
- **AI summaries in the report (optional).** With `generate_ai_summary` and
  `per_section_ai_summaries`, the FGCZ-internal language model writes a summary at
  the top of the MultiQC report and short notes in each section, as reading aids.
- **fastp preprocessing (hidden option).** Setting `max_len1` above 0 through
  `specialOptions` runs the full ezRun fastp preprocessing (length capping, adapter
  and quality trimming) before FastQC. The trimmed reads are used for this QC only.

## How to use it

- **Input dataset:** needs `Name` and `Read1`; `Read2` for paired-end data. SUSHI sets
  `paired` automatically from the presence of `Read2`.
- **Jobs:** one job for the whole dataset.
- **Parameters:**

  | Parameter | Effect |
  |---|---|
  | `paired` | Whether `Read2` files are assessed. SUSHI sets it automatically. |
  | `showNativeReports` | Keeps the per-file FastQC reports and renders the native summary. |
  | `generate_ai_summary` | Add an AI-written summary at the top of the MultiQC report. |
  | `per_section_ai_summaries` | Add AI-written notes to each MultiQC section. |
  | `cmdOptions` | Extra options passed as-is to FastQC. |
  | `specialOptions` | Free-form `key=value` ezRun parameters, e.g. `max_len1` (see above). |
  | `cores`, `ram`, `scratch` | Compute resources. |

- **Outputs:**
  - `MultiQC Report [Link]`: the main report, `multi_FastQC/multiqc_report.html`.
  - `MultiQC [File]`: the MultiQC folder, including its data tables.
  - `FastQC [File]`: per-file FastQC results. The native HTML reports are kept only
    with `showNativeReports`.
  - `FastQC Report [Link]`: only with `showNativeReports`.
- **Gotchas:**
  - Two input files with the same base name produce the same report name and stop
    the run.
  - `perLibrary` appears in the ezRun defaults but has no effect.
  - Log lines such as `Unable to remove path: .../.multiqc_tmp/...` are harmless
    clean-up noise.

## What to describe in the Methods

A Methods description of this step is short, usually two or three sentences.

**Worth describing**

- **FastQC and its version.** FastQC does not print its version; it is the one in the
  job script's `module add` line (`QC/FastQC/<version>`).
- **Single-end or paired-end reads.**
- **What FastQC assessed:** the checks listed under "What it can do".
- **The in-house adapter list** used for adapter content. Name it as such, not by
  its file path.
- **Options passed through `cmdOptions`,** where they change what FastQC assesses.
- **Subsampling, when it happened.** The subsampled input files on the FastQC
  command line end in `-subsample_R1.fastq.gz`.
- **fastp preprocessing, when it happened,** with the options it was given. fastp
  then also applies its built-in quality filter, and may trim poly-G tails on data
  from two-colour Illumina instruments.
- **MultiQC and its version,** used to aggregate the results. The version is on
  MultiQC's banner line (`/// MultiQC ... v<version>`). A `version_check ... now
  available!` line only names a newer release, not the one that ran.

**Not worth describing**

- **FastQC's housekeeping options:** `--extract`, `-o`, `--dir`, `-q`, `-t`
  (threads) and the redirection of its output. They do not change the result.
- **K-mers.** The `--kmers` option has no effect, because FastQC's k-mer module is
  switched off in the FGCZ installation.
- **The AI summaries and the model that writes them.** They are reading aids inside
  the report, not part of the analysis.
- **The native reports,** rendered or not: presentation only.
- **Picard and samtools:** loaded, but never called.
- **Concatenation, read counting and file clean-up.**
- **Counts** of files, reports, report sections or reads. These are results.
- **R, ezRun and SUSHI.** They run the tools; they do not analyse the data.

**Example Methods paragraph** (base case, placeholders in angle brackets):

> Read quality was assessed with FastQC <version> on <single-end or paired-end>
> reads, covering per-base and per-sequence quality, GC and N content, sequence
> length, duplication, overrepresented sequences and adapter content, with an
> in-house list of adapter sequences. The results were summarised with MultiQC
> <version>.
