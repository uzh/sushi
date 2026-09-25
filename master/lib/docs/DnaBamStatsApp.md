# DnaBamStatsApp

Quality control of DNA-based alignments.
ezRun backend: `EzAppDnaBamStats` (`R/app-DnaBamStats.R`).

## Scope

- **Analysis type:** quality control of aligned DNA-based reads (whole-genome or
  targeted DNA sequencing, ChIP-seq, ATAC-seq, CUT&RUN). It gives one report over all
  samples: coverage, mapping quality, duplication, fragment sizes and library
  complexity, with per-sample QC flags.
- **Data:** BAM files from a genome alignment, typically from `Bowtie2App` or
  `BWAApp`.
- **Output:** a report only. The BAM files are not changed.
- **Not covered here:** RNA-seq alignments (`RnaBamStatsApp`); read-level QC
  (`FastqcApp`).
- **Place in a pipeline:** after alignment. No other app takes its output as input.

## What it can do

- **Qualimap `bamqc` per sample (`runQualimap`):** coverage depth and breadth,
  mapping quality, GC content, insert sizes and duplication rate estimates. With
  several samples, a Qualimap `multi-bamqc` summary is also produced.
- **Duplicate metrics with Picard MarkDuplicates (`runPicard`):**
  - duplicate and optical duplicate rates, with the optical duplicate distance set by
    `pixelDist`
  - existing metrics from the aligner (`DupMetrics`, as written by `Bowtie2App`) are
    reused when they were computed with the same distance
- **For paired-end data:**
  - a fragment size distribution, computed with ATACseqQC
  - a library complexity curve, estimated with preseqR from ATACseqQC's
    duplicate-frequency histogram, using a seed fixed in the code
- **Multimapping statistics** from the aligner's tags.
- **Per-sample QC flags:** each sample's metrics are compared with fixed thresholds,
  and the report lists which flags each sample triggers.

## How to use it

- **Input dataset:** needs `Name`, `BAM`, `BAI`, `refBuild`, `Species` and
  `Read Count`. `refFeatureFile` and `paired` are taken from the dataset when
  present.
- **Jobs:** one job for the whole dataset.
- **Parameters:**

  | Parameter | Effect |
  |---|---|
  | `refBuild`, `refFeatureFile` | Reference genome and annotation. |
  | `paired` | Paired-end statistics: fragment sizes and library complexity. |
  | `runQualimap` | Run Qualimap `bamqc`. |
  | `runPicard` | Collect Picard duplicate metrics. |
  | `pixelDist` | Optical duplicate distance for Picard (2500 suits patterned flow cells). |
  | `cores`, `ram`, `scratch` | Compute resources. |

- **Outputs:** `Html [Link]` (the report) and `Report [File]` (the report folder,
  with the Qualimap reports).
- **Gotchas:**
  - Picard metrics from the aligner are reused only when their optical duplicate
    distance equals `pixelDist`.

## What to describe in the Methods

A Methods description of this step is usually one or two sentences. It is a QC step.

**Worth describing**

- **Qualimap and its version** (from the job script's `module add` line), used for
  coverage, mapping quality, GC content and insert sizes. `bamqc` runs per sample,
  and `multi-bamqc` when it ran.
- **Duplicate assessment:** Picard MarkDuplicates metrics with the optical duplicate
  distance (`pixelDist`), whether run here or reused from the alignment step.
- **For paired-end data:** fragment size distributions (ATACseqQC) and library
  complexity extrapolation (preseqR). Their versions appear in the R session listing
  at the end of the log.
- **ezRun and its version,** which assembled the statistics and applied the QC
  thresholds, with R and its version.

**Not worth describing**

- **Qualimap housekeeping:** `-nt` (threads), `-outdir`, `-c`, the `JAVA_OPTS`
  memory setting, and `unset DISPLAY`.
- **The report tables and QC flags,** as presentation of the same metrics.
- **samtools:** a support module.
- **Counts** of reads, coverage values, duplication rates or flags. These are results.
- **SUSHI:** it only launches the job.

**Example Methods paragraph** (paired-end base case, placeholders in angle brackets):

> Alignment quality was assessed with Qualimap <version> (bamqc), covering coverage,
> mapping quality, GC content and insert size. Duplication was assessed from Picard
> <version> MarkDuplicates metrics with an optical duplicate distance of <n> pixels.
> Fragment size distributions were computed with ATACseqQC <version>, and library
> complexity was extrapolated with preseqR <version>, in ezRun <version> under R
> <version>.
