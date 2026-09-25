# MageckTestApp

Gene-level testing of CRISPR screens between two conditions with MAGeCK.
ezRun backend: `EzAppMageckTest` (`R/app-MageckTestApp.R`), report template
`inst/templates/MageckTest.qmd`.

## Scope

- **Analysis type:** comparison of sgRNA counts between a sample group and a
  reference group of a pooled CRISPR screen. The robust rank aggregation test of
  MAGeCK gives positively and negatively selected genes. Optional additions are copy
  number correction, a MAGeCK MLE analysis against a day-0 baseline, and enrichment
  of the hit genes.
- **Data:** per-sample sgRNA count tables from `MageckCountApp`, with a `[Factor]`
  column defining the groups.
- **Output:** gene and sgRNA summary tables, MAGeCK's PDF report, and a report with
  QC, hit lists and enrichment.
- **Not covered here:** counting (`MageckCountApp`, `CountSpacerApp`); exploration of
  counts (`ExploreMageckCountsApp`).
- **Place in a pipeline:** after `MageckCountApp`. It is an end point.

## What it can do

- **Merging.** The count tables of all samples are joined by sgRNA.
- **Copy number and proximity correction (`useCRISPRcleanR`, optional).**
  CRISPRcleanR corrects gene-independent effects of copy number and genomic proximity
  on the counts. The guides are first placed on the genome of `refBuild` with Bowtie
  2. When this correction is applied, the control sgRNAs are not used afterwards.
- **Test with MAGeCK:**
  - `test` of the sample group (`-t`) against the reference group (`-c`), with robust
    rank aggregation at gene level
  - normalisation by `normalizationMethod`: `median` (default), `control` (on the
    control sgRNAs), `total` or `none`
  - gene log fold changes by `geneLFCMethod`
  - the library's control sgRNAs, when `useControls` is on
  - current code runs MAGeCK2 (`mageck2`); earlier runs used MAGeCK (`mageck`)
- **MLE against a day-0 baseline (`day0Label`, optional).** A MAGeCK MLE analysis of
  the conditions relative to day 0. It gives gene effects for both conditions and a
  nine-square plot (MAGeCKFlute).
- **Screen quality.** Separation of the log fold changes of essential and
  non-essential genes (Cohen's d), and recovery of positive control genes when given
  (`positiveControlGenes`).
- **Enrichment (human and mouse):**
  - with `runEnrichment`: over-representation of the hit genes (below `fdrThreshold`)
    in GO biological process (clusterProfiler `enrichGO`) and KEGG pathways
    (`enrichKEGG`)
  - with `runGSEA`: gene set enrichment with fgsea on the MSigDB Hallmark (H) and
    curated (C2) collections
  - with `runPathview`: KEGG pathway maps with pathview

## How to use it

- **Input dataset:** the output rows of `MageckCountApp`: needs `Name`, `Count` and
  `libName`, plus a `[Factor]` column for `grouping`.
- **Jobs:** one job per comparison.
- **Parameters:**

  | Parameter | Effect |
  |---|---|
  | `grouping`, `sampleGroup`, `refGroup` | The factor column and the two groups compared. |
  | `normalizationMethod`, `geneLFCMethod` | MAGeCK normalisation and gene log fold-change method. |
  | `useControls` | Use the library's control sgRNAs. |
  | `useCRISPRcleanR`, `refBuild` | Copy number and proximity correction, and the genome for guide placement. |
  | `day0Label` | Day-0 condition for the MLE analysis. |
  | `fdrThreshold`, `nTopGenes` | Hit threshold, and number of top genes shown. |
  | `runEnrichment`, `runGSEA`, `runPathview`, `species` | Enrichment analyses, for `hsa` or `mmu`. |
  | `positiveControlGenes` | Genes expected as hits, for QC. |
  | `cmdOptions` | Extra options passed as-is to MAGeCK. |
  | `cores`, `ram`, `scratch` | Compute resources. |

- **Outputs:** the gene and sgRNA summaries, the merged count table, MAGeCK's PDF
  report, and the report (`Static Report`).

## What to describe in the Methods

MAGeCK's command runs without an ezRun command-line entry. MAGeCK prints its version
("Welcome to MAGeCK v<version>" or "Welcome to MAGeCK2 v<version>") and its full
command ("Parameters: ... test ...") in the log. That line is the one to follow for
the tool and its version. The job script's `param[['...']]` lines hold the other
settings.

**Worth describing**

- **MAGeCK (or MAGeCK2) and its version,** with the test: `sampleGroup` against
  `refGroup`, robust rank aggregation at gene level.
- **Normalisation** (`median`, `control`, `total` or `none`) and the gene log
  fold-change method, and the use of the control sgRNAs (`--control-sgrna` on the
  command line).
- **CRISPRcleanR correction, when it was used,** with Bowtie 2 placing the guides on
  the genome.
- **The MLE analysis against day 0, when `day0Label` was set.**
- **The hit threshold** (`fdrThreshold`).
- **Enrichment, when it ran:** clusterProfiler for GO biological process and KEGG,
  fgsea on MSigDB Hallmark and C2, and pathview; with their versions from the R
  session listing.
- **Screen quality:** separation of essential and non-essential genes (Cohen's d).

**Not worth describing**

- **The merging of the count tables and `-n` / `--pdf-report`.**
- **The nine-square and other plots:** presentation of the same results.
- **Bowtie 2 when CRISPRcleanR was not used:** loaded, but not called.
- **Counts** of hits, genes or enriched terms. These are results.
- **R, ezRun and SUSHI.** They run the tools; the statistics are MAGeCK's and the
  enrichment packages'.

**Example Methods paragraph** (base case, placeholders in angle brackets):

> Gene-level selection between <sampleGroup> and <refGroup> was tested with MAGeCK
> <version> (mageck test, robust rank aggregation), with <normalisation> normalisation
> <using the library's control sgRNAs> and <gene LFC method> gene log fold changes.
> Genes with a false discovery rate below <fdrThreshold> were considered hits and were
> tested for enrichment of GO biological process terms and KEGG pathways with
> clusterProfiler <version>; gene set enrichment on MSigDB Hallmark and C2 gene sets
> was performed with fgsea <version>.
