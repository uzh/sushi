# CountSpacerApp

Counting of sgRNA spacers in CRISPR screen reads.
ezRun backend: `EzAppCountSpacer` (`R/app-countSpacer.R`), report template
`inst/templates/CountSpacer.qmd`.

## Scope

- **Analysis type:** extraction of the sgRNA spacer sequence from each read of a
  CRISPR screen library, and counting of reads per sgRNA of the screening library,
  per sample.
- **Data:** single-end (Read1) FASTQ reads of pooled CRISPR screens (knock-out,
  CRISPRa, CRISPRi), with a sgRNA library available in the FGCZ library collection.
- **Output:** per-sample sgRNA counts, count statistics and a report.
- **Not covered here:**
  - QC across the samples of a screen (`CrisprScreenQCApp`)
  - MAGeCK counting (`MageckCountApp`) and testing (`MageckTestApp`)
- **Place in a pipeline:** after demultiplexing. Its counts feed `CrisprScreenQCApp`.

## What it can do

- **Read preprocessing with fastp,** with adapter trimming and the usual trimming and
  filtering options. By default it applies an average-quality filter and a minimum
  read length.
- **Spacer extraction:**
  - Reads are kept when the flanking sequences left and right of the spacer are
    found, with up to `maxMismatch` mismatches.
  - The flanking patterns are given (`leftPattern`, `rightPattern`) or, with
    `guessPatterns`, guessed from the reads.
  - The spacer length is given (`spacerLength`), or taken from the library when all
    its sgRNAs have the same length.
  - Extracted spacers shorter than a minimum length are dropped.
- **Alignment with Bowtie (version 1)** of the extracted spacers against the Bowtie
  index of the chosen sgRNA library (`dictPath`), with Bowtie's built-in alignment
  settings.
- **Counting:** reads per sgRNA, with the number of mismatches of each alignment
  (0 to 3) summarised, and counts per target gene.

## How to use it

- **Input dataset:** needs `Name` and `Read1`.
- **Jobs:** one job per sample.
- **Parameters:**

  | Parameter | Effect |
  |---|---|
  | `dictPath` | The sgRNA library from the FGCZ collection. |
  | `leftPattern`, `rightPattern` | Flanking sequences of the spacer; empty lets them be guessed. |
  | `guessPatterns` | Guess the flanking patterns from the reads. |
  | `spacerLength` | Spacer length; 0 takes it from the library. |
  | `maxMismatch` | Mismatches allowed in each flanking pattern. |
  | `trimAdapter` and the fastp parameters | Read preprocessing. A value of 0 switches an option off. |
  | `cores`, `ram`, `scratch` | Compute resources. |

- **Outputs (per sample):** `Count [Link]` (`-sgRNA_counts.txt`), `Html [Link]` (the
  report) and `Report [File]`, with the result table per sgRNA and per target gene.
- **Gotchas:**
  - The flanking patterns actually used, including guessed ones, are shown in the
    report, not in the job log.
  - A run stops when no read aligns to the library, which usually means wrong
    patterns or the wrong library.

## What to describe in the Methods

The job script's `param[['...']]` lines hold the settings. fastp and Bowtie command
lines are logged as `EXECUTED CMD: ...`, and their versions are in the job script's
`module add` line.

**Worth describing**

- **fastp preprocessing,** with its version and the options it was given that are
  switched on, including adapter trimming. fastp also applies its built-in quality
  filter, and may trim poly-G tails on data from two-colour Illumina instruments.
- **Spacer extraction:** by the flanking sequences (given or guessed), with the
  number of mismatches allowed (`maxMismatch`), and the spacer length.
- **Bowtie and its version,** used to align the spacers to the sgRNA library, with its
  built-in alignment settings.
- **The sgRNA library** (`dictPath`), by its name.
- **Counting** of reads per sgRNA and per target gene.

**Not worth describing**

- **Bowtie housekeeping:** `-f` (FASTA input), `-p` (threads), and the output
  processing (`cut`, `sort`).
- **fastp options set to 0,** which are switched off, and fastp's `--thread` and
  `--compression`.
- **The read-count statistics and mismatch summaries.** These are results.
- **R, ezRun and SUSHI.** They run the tools; the counting itself is simple
  tabulation.

**Example Methods paragraph** (placeholders in angle brackets):

> Reads were preprocessed with fastp <version> (<switched-on options>). sgRNA spacers
> were extracted between the flanking sequences <left> and <right>, allowing <n>
> mismatches per flank, aligned with Bowtie <version> to the <library> sgRNA library,
> and counted per sgRNA and per target gene.
