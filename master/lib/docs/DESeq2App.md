# DESeq2App

Differential expression between two groups with DESeq2.
ezRun backend: `EzAppDeseq2` (`R/app-deseq2TwoGroups.R`), with the comparison logic
in `R/twoGroups.R` and the report template `inst/templates/twoGroups.Rmd`.

## Scope

- **Analysis type:** differential expression of genes (or transcripts) between one
  sample group and one reference group, using DESeq2's negative binomial model. An
  optional second factor can be added as a covariate. The results are followed by
  Gene Ontology enrichment.
- **Data:** count tables from `FeatureCountsApp` or `KallistoApp`, with a `[Factor]`
  column defining the groups.
- **Output:** a report with the result table, plots and enrichment results, and a
  data file for the exploreDE Shiny app.
- **Not covered here:** other test frameworks (`EdgeRApp`, `LimmaApp`); exploration
  without testing (`CountQCApp`).
- **Place in a pipeline:** after counting and `CountQCApp`. It is an end point; no
  other app takes its output as input.

## What it can do

All analysis runs in R, in ezRun's code calling DESeq2 and clusterProfiler.

- **Sample selection.** The samples of the two groups being compared are used. With
  `onlyCompGroupsHeatmap` off, the other samples still appear in the heatmaps.
- **Gene filtering.** A gene is tested when it is present (count above a signal
  threshold) in at least half of the samples of either group.
- **Normalisation.** DESeq2 median-ratio size factors, estimated from the present
  genes only.
- **Model and test:**
  - The design is `~ grouping`, or `~ grouping + grouping2` when a second factor is
    given.
  - DESeq2's `DESeq()` fits the negative binomial model with dispersion estimation
    and tests with the Wald test.
  - Outlier replacement is off, and the Cook's distance cutoff is off.
  - Log2 fold changes can be shrunk with ashr (`useLfcShrink`, not on the form).
- **Multiple-testing correction.** False discovery rates are computed by ezRun, with
  Benjamini-Hochberg (`p.adjust`, method `fdr`), over the tested genes only.
- **Enrichment (`runGO`):**
  - Over-representation analysis with clusterProfiler `enricher`, on the GO
    biological process, molecular function and cellular component annotation. The
    candidate genes pass the p-value cut-off (`pValThreshGO`) and the log2 fold-change
    cut-off (`log2RatioThreshGO`). Terms are kept below the adjusted p-value cut-off
    `fdrThreshORA`.
  - Gene set enrichment analysis with clusterProfiler `GSEA`, on the present genes
    ranked by `rankMetric`. Terms are kept below `fdrThreshGSEA`.
  - For several mammalian species at gene level, the up- and down-regulated genes
    (same cut-offs as ORA) are also sent to the Enrichr web service. The top results
    per Enrichr library are shown in the report, next to links to the Enrichr website.
- **Optional RUV correction (`runRUV`, not on the form).** RUVSeq `RUVs` estimates
  `kRUVFactors` factors of unwanted variation from the replicate groups, and the
  corrected counts replace the raw counts before the analysis.
- **Report:** volcano plots, heatmaps of significant genes, and the result table.

## How to use it

- **Input dataset:** needs `Name`, `Count`, `Species`, `refBuild`, `featureLevel` and
  `refFeatureFile`, plus a `[Factor]` column for `grouping`.
- **Jobs:** one job per comparison.
- **Parameters:**

  | Parameter | Effect |
  |---|---|
  | `grouping` | The factor column that defines the groups. |
  | `sampleGroup`, `refGroup` | The two groups compared; fold changes are sample over reference. |
  | `grouping2` | Optional second factor, added to the model as a covariate. |
  | `featureLevel` | `gene` or `isoform`, matching the counts. |
  | `transcriptTypes` | Feature types kept. |
  | `runGO` | Run ORA and GSEA. |
  | `pValThreshGO`, `log2RatioThreshGO` | Cut-offs for ORA candidate genes. |
  | `fdrThreshORA`, `fdrThreshGSEA` | Adjusted p-value cut-offs for enriched terms. |
  | `rankMetric` | Ranking for GSEA: `log2Ratio`, `pValue` or `signedPValue`. |
  | `backgroundExpression` | Offset added before log2 in heatmaps. |
  | `onlyCompGroupsHeatmap` | Show only the compared samples in heatmaps. |
  | `Rversion` | R version loaded for the job. |
  | `cores`, `ram`, `scratch` | Compute resources. |

- **Outputs:** `Static Report [Link]` (the report), `Live Report [Link]` (exploreDE)
  and `Report [File]`.
- **Gotchas:**
  - `sampleGroup` and `refGroup` must differ.
  - The FDR column in the result table is ezRun's, computed over the tested genes. It
    is not DESeq2's own adjusted p-value.

## What to describe in the Methods

The analysis settings are the job script's `param[['...']]` lines. DESeq2's own
progress messages ("using pre-existing size factors", "estimating dispersions",
"fitting model and testing") appear in the log. Package versions appear in the R session
listing at the end of the log.

**Worth describing**

- **DESeq2 and its version, in R and its version.**
- **The comparison:** `sampleGroup` against `refGroup` in the factor `grouping`, and
  the design formula, including `grouping2` as a covariate when present.
- **Gene filtering:** genes tested when present in at least half of the samples of
  either group.
- **Normalisation** by DESeq2 median-ratio size factors estimated from those genes.
- **The test:** negative binomial model with the Wald test, without outlier
  replacement or Cook's distance filtering.
- **Log2 fold-change shrinkage with ashr, when it was used** (the log has a
  `use lfcShrink` line).
- **Multiple-testing correction:** Benjamini-Hochberg false discovery rate over the
  tested genes.
- **Enrichment, when `runGO` is true:**
  - clusterProfiler and its version, for ORA and GSEA on GO BP, MF and CC
  - the candidate-gene cut-offs for ORA (`pValThreshGO`, `log2RatioThreshGO`)
  - the ranking for GSEA (`rankMetric`)
  - the adjusted p-value cut-offs for terms (`fdrThreshORA`, `fdrThreshGSEA`)
  - Enrichr, where it was queried
- **RUV correction, when it was used.**
- **The annotation,** by source and release, from `refBuild`, which follows
  `<organism>/<source>/<genome build>/Annotation/Release<_><release>-<date>`. The
  `-<date>` is when FGCZ set the reference up; it is not part of the release.
- **The feature types kept** (`transcriptTypes`), when restricted.

**Not worth describing**

- **The report's plots and tables, the exploreDE link and the `.h5ad` export:**
  presentation of the same results.
- **`backgroundExpression`:** it only affects the heatmaps.
- **samtools:** a support module.
- **Counts** of significant genes or enriched terms. These are results.
- **SUSHI:** it only launches the job. ezRun's own role (the filtering and the FDR)
  is described through the points above.

**Example Methods paragraph** (two-factor case with GO, placeholders in angle
brackets):

> Differential expression between <sampleGroup> and <refGroup> was tested with DESeq2
> <version> in R <version>, using the design ~ <grouping> + <grouping2>. Genes
> present in at least half of the samples of either group were tested; size factors
> were estimated by the median-ratio method from these genes, and significance was
> assessed with the Wald test, without outlier replacement or Cook's distance
> filtering. P-values were adjusted with the Benjamini-Hochberg method. Gene Ontology
> enrichment was assessed with clusterProfiler <version>, by over-representation
> analysis of genes with p < <pValThreshGO> and by gene set enrichment analysis of
> genes ranked by <rankMetric>, keeping terms with an adjusted p-value below
> <threshold>.
