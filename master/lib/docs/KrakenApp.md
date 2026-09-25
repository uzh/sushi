# KrakenApp

Taxonomic classification of sequencing reads with Kraken 2.
ezRun backend: `EzAppKraken` (`R/app-kraken.R`).

## Scope

- **Analysis type:** assignment of each read (or read pair) to a taxon with Kraken 2,
  giving a taxonomic profile per sample.
- **Data:** FASTQ reads, single-end or paired-end: metagenomic or metatranscriptomic
  samples, amplicon data (with the 16S databases), or host-depleted reads.
- **Output:** a Kraken 2 report and an interactive Krona chart per sample, plus a
  link to the exploreMetaTax Shiny app. Optionally the per-read assignments and the
  unclassified reads.
- **Not covered here:** abundance re-estimation from Kraken reports (`BrackenApp`);
  a quick contamination check of any library, which `FastqScreenApp` already
  includes.
- **Place in a pipeline:** after demultiplexing, often after host-read removal. Its
  reports feed `BrackenApp`.

## What it can do

- **Read preprocessing with fastp** before classification: adapter trimming
  (switchable), front and tail trimming, sliding-window quality trimming, an
  average-quality filter, length capping, poly-X trimming and a minimum read length.
  The trimmed reads are used only for classification and are not kept.
- **Classification with Kraken 2** (`k2 classify`) against one of the databases in
  the FGCZ collection:
  - NCBI-based: Standard, Standard-8, PlusPF, Viral, core_nt, MiniKraken
  - GTDB genome representatives
  - human gut collections: UHGG, HRGM, HROM
  - EuPathDB
  - 16S databases: SILVA 138, Greengenes, RDP
- **Several databases in one pass.** With `multiDB`, several NCBI-taxonomy databases
  can be combined in one classification.
- **Classification settings:** confidence threshold, minimum base quality, and
  minimum number of hit groups.
- **Report extras:** optional per-clade minimizer columns in the report (for
  filtering in exploreMetaTax), optional saving of unclassified reads, and optional
  saving of per-read assignments.
- **Krona chart** per sample, built from the Kraken 2 report.
- **Exclusive mode.** The whole dataset runs on one reserved node, and the database
  is loaded into memory once and shared across samples. The results are the same as
  one job per sample.

## How to use it

- **Input dataset:** needs `Name` and `Read1`; `Read2` when `paired` is on.
- **Jobs:** one job per sample by default; one job for the dataset with `exclusive`.
- **Memory:** a job needs a little more RAM than the size of the chosen database, which
  the form shows next to each database name.
- **Parameters:**

  | Parameter | Effect |
  |---|---|
  | `krakenDBOpt` | Kraken 2 database; several with `multiDB`. |
  | `multiDB` | Classify against all selected databases at once (NCBI-taxonomy databases only). Off: only the first is used. |
  | `paired` | Classify read pairs together. |
  | `krakenConfidenceOpt` | Confidence score threshold, 0 to 1. |
  | `krakenPhredOpt` | Minimum base quality used in classification. |
  | `minimum_hit_groups` | Minimum hit groups needed to make a call. |
  | `report_minimizer_data` | Add minimizer columns to the report. Switched off automatically with several databases. |
  | `save_unclassified`, `save_read_assignments` | Whether the unclassified reads and the per-read assignments are kept. |
  | `trimAdapter` and the fastp parameters | Read preprocessing. A value of 0 switches an option off. |
  | `cmdOptions` | Extra options passed as-is to `k2 classify`. |
  | `exclusive` | One reserved node for the whole dataset, database loaded once. |
  | `cores`, `ram`, `scratch` | Compute resources. |

- **Outputs:** per sample `KrakenReport [File]` (`<sample>.report.txt`),
  `KronaReport [Link]` (`<sample>.html`) and `Live Report [Link]` (exploreMetaTax).
  With the options on, also `KrakenAssignments [File]` and the unclassified reads.
- **Gotchas:**
  - Multi-database classification works only for databases sharing the NCBI
    taxonomy. GTDB, SILVA, UHGG and similar cannot be mixed with them.
  - Several database folder names carry a release date (for example
    `k2_core_nt_20251015`), but many (`Standard`, `PlusPF`, `UHGG`, ...) do not.

## What to describe in the Methods

**Worth describing**

- **fastp preprocessing,** with its version and the options it was given that are
  switched on, including whether adapter trimming was applied. On the fastp command
  line, `--adapter_fasta` means adapter trimming was on and
  `--disable_adapter_trimming` means it was off. fastp also applies its built-in
  quality filter, and may trim poly-G tails on data from two-colour Illumina
  instruments.
- **Kraken 2 and its version,** with its classification settings in words:
  - the confidence threshold (`--confidence`)
  - the minimum number of hit groups (`--minimum-hit-groups`)
  - the minimum base quality (`--minimum-base-quality`)
  - whether read pairs were classified together (`--paired`)
  - any options added through `cmdOptions`
- **The database,** by name, not by path. It is the folder on the `--db` path:
  - A release date is part of some folder names (`k2_core_nt_20251015` is core_nt of
    2025-10-15). Many folder names carry none (`Standard`, `PlusPF`, `UHGG`, ...).
  - Several comma-separated folders mean several databases were combined in one
    classification.
- **Krona and its version,** used to draw an interactive chart of each sample's
  classification. A short mention is enough.

**Not worth describing**

- **Housekeeping options of `k2 classify`:** `--use-daemon`, `--use-names`,
  `--output`, `--report`, `--unclassified-out`, `--threads`, the log redirection, and
  the `k2 clean --stop-daemon` calls.
- **`--report-minimizer-data`:** it only adds columns to the report.
- **Exclusive mode and the database daemon.** They change how the job ran, not the
  result.
- **fastp options set to 0,** which are switched off, and fastp's `--thread` and
  `--compression`.
- **Compression (`pigz`), deletion of the trimmed reads, and the conversion of the
  report into Krona's input format.**
- **Counts** of reads, classified reads or taxa. These are results.
- **R, ezRun and SUSHI.** They run the tools; they do not analyse the data.

**Example Methods paragraph** (placeholders in angle brackets):

> Reads were preprocessed with fastp <version> (<switched-on options>) and
> classified <as read pairs> with Kraken 2 <version> against the <database> database
> <of <release date>>, with a confidence threshold of <confidence>, a minimum of <n>
> hit groups and a minimum base quality of <q>. Classification results were
> visualised with Krona <version>.
