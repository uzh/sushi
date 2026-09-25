# MageckCountApp

Counting of sgRNAs in CRISPR screen reads with MAGeCK2.
ezRun backend: `EzAppMageckCount` (`R/app-MageckCountApp.R`), report template
`inst/templates/MageckCountQC.qmd`.

## Scope

- **Analysis type:** counting of reads per sgRNA of a screening library, per sample,
  with `mageck2 count`, including a count QC report.
- **Data:** FASTQ reads (Read1) of pooled CRISPR screens, with a sgRNA library
  available in the FGCZ library collection.
- **Output:** per-sample sgRNA count table, count summary and QC report.
- **Not covered here:** counting by flanking-pattern extraction and Bowtie
  (`CountSpacerApp`); library identification (`CrisprScreenQCApp`); testing
  (`MageckTestApp`).
- **Place in a pipeline:** after demultiplexing. Its count tables feed
  `ExploreMageckCountsApp` and `MageckTestApp`.

## What it can do

- **Library files.** The chosen library (`libName`) provides the MAGeCK sgRNA table
  and, when available, a list of control sgRNAs. When the library has no MAGeCK files
  yet, they are generated from the library table (sgRNA, sequence, gene, control
  flag).
- **Counting with `mageck2 count`:**
  - reads from the raw FASTQ files, without separate trimming
  - MAGeCK2 determines the 5' trim length of the reads itself
  - the control sgRNAs are passed as `--control-sgrna` when the library has them
  - further options can be added through `cmdOptions`
- **QC report:** from MAGeCK2's count summary (mapped reads, zero-count sgRNAs, Gini
  index and similar).

## How to use it

- **Input dataset:** needs `Name` and `Read1`.
- **Jobs:** one job per sample.
- **Parameters:**

  | Parameter | Effect |
  |---|---|
  | `libName` | The sgRNA library from the FGCZ collection. |
  | `cmdOptions` | Extra options passed as-is to `mageck2 count`. |
  | `cores`, `ram`, `scratch` | Compute resources. |

- **Outputs (per sample):** the count table (`<sample>.count.txt`), the count summary
  (`<sample>.countsummary.txt`) and the QC report (`<sample>.html`).

## What to describe in the Methods

The `mageck2 count` command line is logged as `EXECUTED CMD: ...`. MAGeCK2 also
prints its version ("Welcome to MAGeCK2 v<version>") and its full command
("Parameters: ...").

**Worth describing**

- **MAGeCK2 and its version** (`mageck2 count`), with automatic determination of the
  5' trim length.
- **The sgRNA library** by name (from `libName`), and the use of its control sgRNAs
  when `--control-sgrna` is on the command line.
- **Options from `cmdOptions`,** in words.

**Not worth describing**

- **The generation of the MAGeCK library files and the QC report.**
- **`-n`** (the output name) and the conda environment activation.
- **Counts** of reads, mapped reads, zero-count sgRNAs or the Gini index. These are
  results.
- **R, ezRun and SUSHI.** They run MAGeCK2; they do not analyse the data.

**Example Methods paragraph** (placeholders in angle brackets):

> Reads were counted per sgRNA with MAGeCK2 <version> (mageck2 count) against the
> <library> sgRNA library, including its control sgRNAs, with the 5' trim length
> determined automatically.
