# EdgeRApp

Differential expression between two groups with edgeR.
ezRun backend: `EzAppEdger` (`R/app-edgerTwoGroups.R`), with the comparison logic in
`R/twoGroups.R` and the report template `inst/templates/twoGroups.Rmd`.

## Scope

- **Analysis type:** differential expression of genes (or transcripts) between one
  sample group and one reference group with edgeR. An optional second factor can be
  added as a covariate, and a difference of differences can be tested against two
  baseline groups. The results are followed by Gene Ontology enrichment.
- **Data:** count tables from `FeatureCountsApp` or `KallistoApp`, with a `[Factor]`
  column defining the groups.
- **Output:** a report with the result table, plots and enrichment results, and a
  data file for the exploreDE Shiny app.
- **Not covered here:** DESeq2 (`DESeq2App`); limma-voom (`LimmaApp`); exploration
  without testing (`CountQCApp`).
- **Place in a pipeline:** after counting and `CountQCApp`. It is an end point; no
  other app takes its output as input.

## What it can do

All analysis runs in R, in ezRun's code calling edgeR and clusterProfiler.

- **Gene filtering.** A gene is tested when it is present (count above a signal
  threshold) in at least half of the samples of either group.
- **Normalisation** with edgeR `calcNormFactors`: `TMM` (default), `RLE`,
  `upperquartile` or `none`.
- **Two test frameworks (`testMethod`):**
  - `glm` (default):
    - design `~ 0 + group`, plus `grouping2` as a covariate when given
    - dispersions estimated with `estimateDisp`
    - with `deTest = QL` (default), the quasi-likelihood F-test (`glmQLFit` and
      `glmQLFTest`); with `LR`, the likelihood-ratio test (`glmFit` and `glmLRT`)
    - when `sampleGroupBaseline` and `refGroupBaseline` are both given, the tested
      contrast is (sample - sample baseline) - (reference - reference baseline)
  - `exactTest`: edgeR's exact test for two groups, after `estimateDisp`. It uses
    no second factor.
- **Log2 fold changes** are computed with a prior count equal to
  `backgroundExpression` (default 10), which shrinks the fold changes of genes with
  low counts.
- **Multiple-testing correction.** False discovery rates are computed by ezRun, with
  Benjamini-Hochberg (`p.adjust`, method `fdr`), over the tested genes only.
- **Enrichment (`runGO`):**
  - Over-representation analysis with clusterProfiler `enricher` on GO BP, MF and CC,
    for candidate genes passing `pValThreshGO` and `log2RatioThreshGO`, with terms
    kept below `fdrThreshORA`.
  - Gene set enrichment analysis with clusterProfiler `GSEA` on the present genes
    ranked by `rankMetric`, with terms kept below `fdrThreshGSEA`.
  - Links to the Enrichr website for the up- and down-regulated genes. Enrichr
    itself is not queried.
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
  | `sampleGroupBaseline`, `refGroupBaseline` | Optional baseline groups for a difference-of-differences contrast (glm only). |
  | `grouping2` | Optional second factor, added to the model as a covariate (glm only). |
  | `testMethod` | `glm` or `exactTest`. |
  | `deTest` | `QL` (quasi-likelihood F-test) or `LR` (likelihood-ratio test), for glm. |
  | `normMethod` | `TMM`, `RLE`, `upperquartile` or `none`. |
  | `backgroundExpression` | Prior count for the log2 fold changes, and the offset in heatmaps. |
  | `featureLevel`, `transcriptTypes` | Feature level, and the feature types kept. |
  | `runGO`, `pValThreshGO`, `log2RatioThreshGO`, `fdrThreshORA`, `rankMetric`, `fdrThreshGSEA` | Enrichment analysis, as above. |
  | `pValueHighlightThresh`, `log2RatioHighlightThresh` | Highlighting thresholds in the plots only. |
  | `Rversion` | R version loaded for the job. |
  | `cores`, `ram`, `scratch` | Compute resources. |

- **Outputs:** `Static Report [Link]` (the report), `Live Report [Link]` (exploreDE)
  and `Report [File]`.
- **Gotchas:**
  - The FDR column in the result table is ezRun's, computed over the tested genes.

## What to describe in the Methods

The analysis settings are the job script's `param[['...']]` lines. Package versions
appear in the R session listing at the end of the log.

**Worth describing**

- **edgeR and its version, in R and its version.**
- **The comparison:** `sampleGroup` against `refGroup` in the factor `grouping`,
  the design, including `grouping2` as a covariate when present, and the baseline
  groups of a difference-of-differences contrast when given.
- **Gene filtering:** genes tested when present in at least half of the samples of
  either group.
- **Normalisation** with `normMethod` (for example "TMM normalisation").
- **The test:** for `glm`, dispersion estimation and the quasi-likelihood F-test or
  the likelihood-ratio test (`deTest`); for `exactTest`, edgeR's exact test.
- **The prior count for the log2 fold changes** (`backgroundExpression`).
- **Multiple-testing correction:** Benjamini-Hochberg false discovery rate over the
  tested genes.
- **Enrichment, when `runGO` is true:** clusterProfiler and its version, ORA and GSEA
  on GO BP, MF and CC, with the candidate-gene cut-offs, the GSEA ranking, and the
  adjusted p-value cut-offs for terms.
- **RUV correction, when it was used.**
- **The annotation,** by source and release, from `refBuild`, which follows
  `<organism>/<source>/<genome build>/Annotation/Release<_><release>-<date>`. The
  `-<date>` is when FGCZ set the reference up; it is not part of the release.
- **The feature types kept** (`transcriptTypes`), when restricted.

**Not worth describing**

- **The Enrichr links:** they only open the website.
- **`pValueHighlightThresh` and `log2RatioHighlightThresh`:** plot highlighting only.
- **The report's plots and tables, the exploreDE link and the `.h5ad` export.**
- **samtools:** a support module.
- **Counts** of significant genes or enriched terms. These are results.
- **SUSHI:** it only launches the job.

**Example Methods paragraph** (glm with QL, placeholders in angle brackets):

> Differential expression between <sampleGroup> and <refGroup> was tested with edgeR
> <version> in R <version>. Genes present in at least half of the samples of either
> group were tested. Counts were normalised with the TMM method, dispersions were
> estimated, and a generalised linear model (design ~ 0 + <grouping>) was tested
> with the quasi-likelihood F-test, computing log2 fold changes with a prior count of
> <backgroundExpression>. P-values were adjusted with the Benjamini-Hochberg method.
> Gene Ontology enrichment was assessed with clusterProfiler <version> by
> over-representation and gene set enrichment analysis.
