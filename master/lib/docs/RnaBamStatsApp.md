# RnaBamStatsApp

Quality control of RNA-seq alignments.
ezRun backend: `EzAppRnaBamStats` (`R/app-RnaBamStats.R`), report template
`inst/templates/RNABamStats.Rmd`.

## Scope

- **Analysis type:** quality control of aligned RNA-seq reads. It gives one report
  over all samples, on alignment, read distribution, coverage, splicing, errors and
  duplication.
- **Data:** BAM files from a genome alignment of bulk RNA-seq reads, usually from
  `STARApp`.
- **Output:** a report only. The BAM files are not changed.
- **Not covered here:** DNA alignments (`DnaBamStatsApp`); QC of count tables
  (`CountQCApp`).
- **Place in a pipeline:** after `STARApp`, next to `FeatureCountsApp`. No other app
  takes its output as input.

## What it can do

The statistics are computed by ezRun's own R code, using the Bioconductor packages
GenomicAlignments and Rsamtools, with these components:

- **Alignment statistics:** how many reads align to one or to several locations
  (from the aligner's multimapping tag), and reads per chromosome.
- **Read distribution:** reads on exons and introns of protein-coding transcripts, on
  the exons of the other transcript types, on RNA repeats (by class and family, where
  the reference provides a repeat annotation), and reads outside any annotation.
- **Gene body coverage:** coverage along transcripts, from a sample of transcripts,
  grouped by expression level.
- **Fragment size distribution** for paired-end data.
- **Position-specific error rates:** mismatch rates along the read, from a sample of
  reads compared with the reference genome.
- **Splice junctions:** annotated versus novel junctions, and junction saturation.
  They are taken from the aligner's junction file when present, and otherwise
  computed from the BAM files.
- **Duplication rate against expression level,** computed with dupRadar. The rates
  computed by `STARApp` are reused when the annotation, strandedness and pairing
  match; otherwise they are recomputed.
- **Strandedness,** read from the RSeQC result produced by `STARApp`, when present.
- **An IGV session link.**

## How to use it

- **Input dataset:** needs `Name`, `BAM`, `BAI`, `refBuild` and `Species`.
  `refFeatureFile`, `paired` and `strandMode` are taken from the dataset when
  present.
- **Jobs:** one job for the whole dataset.
- **Parameters:**

  | Parameter | Effect |
  |---|---|
  | `refBuild`, `refFeatureFile` | Reference genome and annotation for the read distribution, coverage and junctions. |
  | `strandMode` | `both` (unstranded), `sense` or `antisense`. |
  | `paired` | Paired-end statistics, including fragment sizes. |
  | `cores`, `ram`, `scratch` | Compute resources. |

- **Outputs:** `Html [Link]` (the report, `00index.html`) and `Report [File]` (the
  report folder).
- **Gotchas:**
  - The upstream duplication rates are reused only when the dataset's
    `refFeatureFile`, `strandMode` and `paired` equal this run's settings. The log
    says which applied: "using the duplication rates computed by the aligner" or
    "recomputing".

## What to describe in the Methods

A Methods description of this step is usually one or two sentences. It is a QC step.

**Worth describing**

- **What was assessed:** alignment and multimapping statistics, read distribution over
  exons, introns, transcript types and repeats, gene body coverage, fragment sizes (paired-end), position-specific
  error rates, splice junctions, duplication rates and strandedness.
- **The software.** ezRun and its version, which computed the statistics in R, with R
  and its version. dupRadar for the duplication rates. The versions of ezRun,
  dupRadar, GenomicAlignments and Rsamtools appear in the R session listing at the
  end of the log.
- **The annotation** used for read distribution and coverage. The `refBuild` value
  follows `<organism>/<source>/<genome build>/Annotation/Release<_><release>-<date>`.
  The `-<date>` is when FGCZ set the reference up; it is not part of the release.
- **Reuse of upstream results,** where the log says so: splice junctions from the
  aligner's junction files, duplication rates from the alignment step, and
  strandedness from RSeQC.

**Not worth describing**

- **The report and the IGV link.**
- **Modules loaded as support:** Picard, BamUtil, samtools and Python are loaded for
  the fallback computations.
- **Counts** of reads, percentages, error rates or duplication rates. These are
  results.
- **SUSHI:** it only launches the job.

**Example Methods paragraph** (placeholders in angle brackets):

> Alignment quality was assessed with ezRun <version> in R <version>, covering
> alignment and multimapping statistics, the distribution of reads over annotated
> feature types of the <source> annotation release <release>, gene body coverage,
> fragment sizes, position-specific error rates and splice junctions. Duplication
> rates relative to expression were computed with dupRadar <version>.
